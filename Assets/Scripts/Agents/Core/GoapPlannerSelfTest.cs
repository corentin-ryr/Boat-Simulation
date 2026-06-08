using System.Collections.Generic;
using Agents.Economy;
using TerrainGrid;
using UnityEngine;

namespace Agents.Core
{
    // Smoke test for the planner — runs three planning scenarios against
    // a fake IGoapAgent + fake action library and Debug.Logs the results.
    // Attach to any GameObject and tick `runOnStart` (or press F7) to
    // invoke. Pure C# tests don't need this — they could ctor the planner
    // directly — but we don't have a Unity-test rig in this project yet,
    // so this is the cheapest standalone validation.
    //
    // Scenarios:
    //   1. Earn flow: Harvest → Bake → Sell. Expected: 3-action plan.
    //   2. Impossible goal: targets a key no action provides. Expected:
    //      null plan (within node cap).
    //   3. Adversarial fan-out: many no-op actions to stress the closed
    //      set. Expected: terminates without exceeding MaxNodes.
    public class GoapPlannerSelfTest : MonoBehaviour
    {
        public bool runOnStart = true;
        public KeyCode runKey = KeyCode.F7;

        void Start()
        {
            if (runOnStart) RunAll();
        }

        void Update()
        {
            if (Input.GetKeyDown(runKey)) RunAll();
        }

        public void RunAll()
        {
            Debug.Log("[GoapPlannerSelfTest] running ...");
            TestEarnFlow();
            TestImpossibleGoal();
            TestNoOpFanOut();
            Debug.Log("[GoapPlannerSelfTest] done.");
        }

        void TestEarnFlow()
        {
            FakeAgent a = new FakeAgent();
            a.State.Set(WorldKey.HasWheat, false);
            a.State.Set(WorldKey.HasBread, false);
            a.State.Set(WorldKey.MoneyEnough, false);

            List<GoapAction> actions = new List<GoapAction>
            {
                new HarvestStub(), new BakeStub(), new SellStub(),
            };
            GoapGoal goal = new EarnStubGoal();

            Queue<GoapAction> plan = GoapPlanner.Plan(a.State, goal, actions, a);
            if (plan == null) { Debug.LogError("Earn flow: plan was null"); return; }
            Debug.Log($"Earn flow: plan length = {plan.Count}, steps = {string.Join(" → ", plan)}");
            if (plan.Count != 3) Debug.LogWarning($"Earn flow: expected 3 steps, got {plan.Count}");
        }

        void TestImpossibleGoal()
        {
            FakeAgent a = new FakeAgent();
            List<GoapAction> actions = new List<GoapAction> { new HarvestStub() };
            GoapGoal goal = new EarnStubGoal();   // requires MoneyEnough — no action provides it here
            Queue<GoapAction> plan = GoapPlanner.Plan(a.State, goal, actions, a);
            Debug.Log($"Impossible goal: plan = {(plan == null ? "null (correct)" : "non-null (bug)")}");
        }

        void TestNoOpFanOut()
        {
            // Builds a 10-action library whose effects all collapse to
            // the same state — closed-set should catch each branch.
            // Goal is MoneyEnough so the planner has to actually run
            // a no-op fanout before giving up at MaxNodes.
            FakeAgent a = new FakeAgent();
            List<GoapAction> actions = new List<GoapAction>();
            for (int i = 0; i < 10; i++) actions.Add(new NoOpStub($"noop{i}"));

            Queue<GoapAction> plan = GoapPlanner.Plan(a.State, new EarnStubGoal(), actions, a);
            Debug.Log($"No-op fan-out: plan = {(plan == null ? "null" : plan.Count.ToString())} (expected null)");
        }

        // -------- Test fixtures --------

        sealed class FakeAgent : IGoapAgent
        {
            public WorldState State { get; } = new WorldState();
            public Vector3 Position { get; set; }
            public Inventory Inventory { get; } = new Inventory();
            public int Money { get; set; }
            public (ChunkCoord Coord, int CellIndex)? AtTile { get; set; }
            public (ChunkCoord Coord, int CellIndex)? Home   => null;
            public (ChunkCoord Coord, int CellIndex)? Job    => null;
            public string DisplayName => "fake";
            readonly Dictionary<System.Type, ActionData> bucket = new Dictionary<System.Type, ActionData>();
            public T GetActionData<T>() where T : ActionData, new()
            {
                System.Type t = typeof(T);
                if (!bucket.TryGetValue(t, out ActionData d)) { d = new T(); bucket[t] = d; }
                return (T)d;
            }
        }

        sealed class HarvestStub : GoapAction
        {
            public override string Name => "Harvest";
            public override IEnumerable<KeyValuePair<string, object>> Preconditions { get; }
                = new List<KeyValuePair<string, object>> { new KeyValuePair<string, object>(WorldKey.HasWheat, false) };
            public override IEnumerable<KeyValuePair<string, object>> Effects { get; }
                = new List<KeyValuePair<string, object>> { new KeyValuePair<string, object>(WorldKey.HasWheat, true) };
            public override float Cost(IGoapAgent _) => 5f;
            public override bool Perform(IGoapAgent _, float dt) => true;
        }

        sealed class BakeStub : GoapAction
        {
            public override string Name => "Bake";
            public override IEnumerable<KeyValuePair<string, object>> Preconditions { get; }
                = new List<KeyValuePair<string, object>>
                {
                    new KeyValuePair<string, object>(WorldKey.HasWheat, true),
                    new KeyValuePair<string, object>(WorldKey.HasBread, false),
                };
            public override IEnumerable<KeyValuePair<string, object>> Effects { get; }
                = new List<KeyValuePair<string, object>>
                {
                    new KeyValuePair<string, object>(WorldKey.HasWheat, false),
                    new KeyValuePair<string, object>(WorldKey.HasBread, true),
                };
            public override float Cost(IGoapAgent _) => 3f;
            public override bool Perform(IGoapAgent _, float dt) => true;
        }

        sealed class SellStub : GoapAction
        {
            public override string Name => "Sell";
            public override IEnumerable<KeyValuePair<string, object>> Preconditions { get; }
                = new List<KeyValuePair<string, object>>
                {
                    new KeyValuePair<string, object>(WorldKey.HasBread, true),
                    new KeyValuePair<string, object>(WorldKey.MoneyEnough, false),
                };
            public override IEnumerable<KeyValuePair<string, object>> Effects { get; }
                = new List<KeyValuePair<string, object>>
                {
                    new KeyValuePair<string, object>(WorldKey.HasBread, false),
                    new KeyValuePair<string, object>(WorldKey.MoneyEnough, true),
                };
            public override float Cost(IGoapAgent _) => 1f;
            public override bool Perform(IGoapAgent _, float dt) => true;
        }

        sealed class NoOpStub : GoapAction
        {
            readonly string name;
            public NoOpStub(string n) { name = n; }
            public override string Name => name;
            public override IEnumerable<KeyValuePair<string, object>> Preconditions { get; }
                = System.Array.Empty<KeyValuePair<string, object>>();
            // Empty effects — closed set should reject the resulting
            // (unchanged) state immediately.
            public override IEnumerable<KeyValuePair<string, object>> Effects { get; }
                = System.Array.Empty<KeyValuePair<string, object>>();
            public override bool Perform(IGoapAgent _, float dt) => true;
        }

        sealed class EarnStubGoal : GoapGoal
        {
            public override string Name => "Earn";
            public override IEnumerable<KeyValuePair<string, object>> TargetState { get; }
                = new List<KeyValuePair<string, object>>
                {
                    new KeyValuePair<string, object>(WorldKey.MoneyEnough, true),
                };
            public override float Priority(IGoapAgent _) => 100f;
        }
    }
}
