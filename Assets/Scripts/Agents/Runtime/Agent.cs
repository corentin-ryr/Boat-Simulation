using System;
using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using TerrainGrid;
using UnityEngine;
using World;

namespace Agents.Runtime
{
    // MonoBehaviour realization of IGoapAgent. Owns per-agent state
    // (Inventory, Money, Needs, Home, Job, plan, action-data buckets)
    // and ticks the GOAP loop on a fixed cadence.
    //
    // Loop structure (split across the wall-clock frame):
    //   Every plan tick (planTickInterval, wall-seconds):
    //     1. Tick sensors (Always / Once / Interval mode-gated).
    //     2. Select goal (highest unsatisfied priority + hysteresis).
    //     3. Ensure plan (re-plan if missing or current action invalid).
    //   Every frame:
    //     4. Execute current action (Perform → done? → Complete + pop).
    //     5. Tick needs (Hunger / Energy decay against in-game hours).
    //
    // Why split: action movement (MoveToTileAction tweening
    // transform.position) needs to read every frame's deltaTime to look
    // smooth, but planning / sensing only need fresh data 4× a second.
    // Decoupling keeps the planner cheap at 100 agents.
    //
    // Initialization contract: spawner instantiates the prefab, calls
    // Initialize(blueprint, name, home, job) once. Until then, Update()
    // is a no-op (Blueprint == null).
    public class Agent : MonoBehaviour, IGoapAgent
    {
        [Header("Identity")]
        [SerializeField] string displayName = "Agent";
        public string DisplayName => displayName;

        [Header("State")]
        public AgentNeeds Needs = new AgentNeeds();
        [SerializeField] int money = 0;
        public int Money { get => money; set => money = value; }
        Inventory inventory = new Inventory();
        public Inventory Inventory => inventory;

        [Header("Tuning")]
        [Tooltip("Wall-clock seconds between plan-tick passes (sensors + goal + plan). " +
                 "Action.Perform still runs every frame so movement stays smooth.")]
        public float planTickInterval = 0.25f;

        [Tooltip("Once a goal is selected, the agent commits to it for at least this many " +
                 "game-time seconds before re-evaluating against newly-bumped priorities. " +
                 "Prevents thrashing between Survive and Earn when both score similarly. " +
                 "Default 300 ≈ 5 game-minutes (10s wall-time at default 30s/game-hour). " +
                 "Set to 0 for pure 'always pick max priority' behaviour.")]
        public float minGoalGameSeconds = 300f;

        // -------- IGoapAgent surface --------
        public WorldState State { get; } = new WorldState();
        public Vector3 Position
        {
            get => transform.position;
            set => transform.position = value;
        }
        public (ChunkCoord Coord, int CellIndex)? AtTile { get; set; }

        (ChunkCoord, int)? home;
        (ChunkCoord, int)? job;
        public (ChunkCoord Coord, int CellIndex)? Home => home;
        public (ChunkCoord Coord, int CellIndex)? Job => job;

        readonly Dictionary<Type, ActionData> actionDataBuckets = new Dictionary<Type, ActionData>();
        public T GetActionData<T>() where T : ActionData, new()
        {
            Type t = typeof(T);
            if (!actionDataBuckets.TryGetValue(t, out ActionData d))
            {
                d = new T();
                actionDataBuckets[t] = d;
            }
            return (T)d;
        }

        // -------- Runtime --------
        public GoapBlueprint Blueprint { get; private set; }
        public GoapGoal CurrentGoal { get; private set; }
        public GoapAction CurrentAction { get; private set; }
        public Queue<GoapAction> CurrentPlan { get; private set; }

        // Sleeping flag is set by SleepAction.Start / cleared on Stop /
        // Complete so AgentNeeds.Tick knows whether to recover Energy
        // instead of decaying it. Not on IGoapAgent because nothing in
        // Core cares — only the agent's own Needs update reads it.
        public bool IsSleeping { get; set; }

        float planTickAccum;
        float goalCommittedGameSeconds;

        // One slot per sensor in Blueprint.Sensors. NextDue is game-time
        // seconds; OnceFired tracks the Once-mode sentinel.
        float[] sensorNextDue;
        bool[] sensorOnceFired;

        public void Initialize(GoapBlueprint blueprint,
                               string newDisplayName,
                               (ChunkCoord Coord, int CellIndex)? newHome,
                               (ChunkCoord Coord, int CellIndex)? newJob)
        {
            Blueprint = blueprint;
            displayName = newDisplayName ?? displayName;
            home = newHome;
            job  = newJob;
            // Allocate sensor scheduling arrays to match the blueprint size.
            int n = Blueprint?.Sensors.Count ?? 0;
            sensorNextDue = new float[n];
            sensorOnceFired = new bool[n];
            goalCommittedGameSeconds = -1f;
        }

        void OnEnable() => AgentRegistry.Instance?.Register(this);
        void OnDisable()
        {
            // Release any reserved building slots — no zombies left behind.
            BuildingRegistry.Instance?.ReleaseAllForOwner(this);
            AgentRegistry.Instance?.Unregister(this);
        }

        void OnDestroy()
        {
            BuildingRegistry.Instance?.ReleaseAllForOwner(this);
            AgentRegistry.Instance?.Unregister(this);
        }

        void Update()
        {
            if (Blueprint == null) return;
            float dt = Time.deltaTime;

            // Needs continuously decay against game-time — F2 fast-forward
            // affects this through GameClock.TimeOfDayHours.
            Needs.Tick(IsSleeping);

            // Plan tick on interval — wall-clock cadence, not game-clock,
            // so 4 Hz of planning regardless of timeScale.
            planTickAccum += dt;
            if (planTickAccum >= planTickInterval)
            {
                planTickAccum = 0f;
                TickSensors();
                SelectGoal();
                EnsurePlan();
            }

            ExecuteAction(dt);
        }

        // -------- Sensor scheduling --------

        void TickSensors()
        {
            if (Blueprint.Sensors.Count == 0) return;
            float gameSeconds = CurrentGameTimeSeconds();
            for (int i = 0; i < Blueprint.Sensors.Count; i++)
            {
                GoapSensor s = Blueprint.Sensors[i];
                switch (s.Mode)
                {
                    case GoapSensor.SenseMode.Always:
                        s.Sense(this);
                        break;
                    case GoapSensor.SenseMode.Once:
                        if (!sensorOnceFired[i])
                        {
                            s.Sense(this);
                            sensorOnceFired[i] = true;
                        }
                        break;
                    case GoapSensor.SenseMode.Interval:
                        if (gameSeconds >= sensorNextDue[i])
                        {
                            s.Sense(this);
                            sensorNextDue[i] = gameSeconds + s.IntervalSeconds;
                        }
                        break;
                }
            }
        }

        // -------- Goal selection --------

        void SelectGoal()
        {
            if (Blueprint.Goals.Count == 0) return;
            float gameSeconds = CurrentGameTimeSeconds();

            // Honour the commit window unless the current goal is already
            // satisfied (no work left) or there's no current goal yet.
            if (CurrentGoal != null && goalCommittedGameSeconds >= 0f)
            {
                bool stillCommitted = (gameSeconds - goalCommittedGameSeconds) < minGoalGameSeconds;
                bool stillUnsatisfied = !State.Satisfies(CurrentGoal.TargetState);
                if (stillCommitted && stillUnsatisfied) return;
            }

            GoapGoal best = null;
            float bestPriority = -1f;
            for (int i = 0; i < Blueprint.Goals.Count; i++)
            {
                GoapGoal g = Blueprint.Goals[i];
                if (State.Satisfies(g.TargetState)) continue;
                float p = g.Priority(this);
                if (p > bestPriority) { bestPriority = p; best = g; }
            }

            if (best != CurrentGoal)
            {
                if (CurrentAction != null) { CurrentAction.Stop(this); CurrentAction = null; }
                CurrentPlan = null;
                CurrentGoal = best;
                goalCommittedGameSeconds = gameSeconds;
            }
        }

        // -------- Plan / execution --------

        void EnsurePlan()
        {
            if (CurrentGoal == null) return;

            bool needReplan = CurrentPlan == null;
            if (CurrentAction != null && !CurrentAction.IsValid(this)) needReplan = true;
            if (!needReplan) return;

            if (CurrentAction != null) { CurrentAction.Stop(this); CurrentAction = null; }
            CurrentPlan = GoapPlanner.Plan(State.Clone(), CurrentGoal, Blueprint.Actions, this);
            if (CurrentPlan != null && CurrentPlan.Count > 0)
            {
                CurrentAction = CurrentPlan.Dequeue();
                CurrentAction.Start(this);
            }
        }

        void ExecuteAction(float dt)
        {
            if (CurrentAction == null) return;

            // Per-frame validity recheck — slot stolen, inventory ran
            // out, target tile became unloaded. Failure invalidates the
            // whole plan; the next plan tick re-plans.
            if (!CurrentAction.IsValid(this))
            {
                CurrentAction.Stop(this);
                CurrentAction = null;
                CurrentPlan = null;
                return;
            }

            bool done = CurrentAction.Perform(this, dt);
            if (!done) return;

            CurrentAction.Complete(this);
            CurrentAction = null;
            if (CurrentPlan != null && CurrentPlan.Count > 0)
            {
                CurrentAction = CurrentPlan.Dequeue();
                CurrentAction.Start(this);
            }
        }

        // -------- Helpers --------

        static float CurrentGameTimeSeconds()
        {
            GameClock c = GameClock.Instance;
            if (c == null) return 0f;
            return (c.DayOfYear * GameClock.HoursPerDay + c.TimeOfDayHours) * 3600f;
        }
    }
}
