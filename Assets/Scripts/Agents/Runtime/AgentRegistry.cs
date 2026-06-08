using System.Collections.Generic;
using TerrainGrid;
using UnityEngine;

namespace Agents.Runtime
{
    // Authoritative list of every live Agent in the scene. MonoBehaviour
    // singleton so other systems (the debug HUD, future quest hooks, AI
    // sensors that need to know neighbours) can iterate without
    // FindObjectsOfType every frame.
    //
    // Agent.OnEnable / OnDisable call Register / Unregister, so the list
    // is a tight reflection of "what's actually playing right now."
    // Despawning an Agent removes it; re-enabling re-registers.
    //
    // At v1 scale (100 agents) linear-scan queries are fine. If we ever
    // push the agent count up enough that "find nearest" becomes a hot
    // path, swap in a chunked spatial hash here without touching callers.
    public class AgentRegistry : MonoBehaviour
    {
        public static AgentRegistry Instance { get; private set; }

        readonly List<Agent> agents = new List<Agent>();
        public IReadOnlyList<Agent> All => agents;
        public int Count => agents.Count;

        void Awake()
        {
            if (Instance != null && Instance != this) { Destroy(this); return; }
            Instance = this;
        }

        void OnDestroy()
        {
            if (Instance == this) Instance = null;
        }

        public void Register(Agent a)
        {
            if (a == null) return;
            if (!agents.Contains(a)) agents.Add(a);
        }

        public void Unregister(Agent a)
        {
            if (a == null) return;
            agents.Remove(a);
        }

        // Linear scan finds (in v1) — cheap at 100 agents, faster than a
        // dictionary because we touch every agent once anyway.
        public Agent FindNearest(Vector3 pos)
        {
            Agent best = null;
            float bestSq = float.PositiveInfinity;
            for (int i = 0; i < agents.Count; i++)
            {
                Agent a = agents[i];
                if (a == null) continue;
                float d = (a.transform.position - pos).sqrMagnitude;
                if (d < bestSq) { bestSq = d; best = a; }
            }
            return best;
        }

        // Agents whose Job tile is the given (coord, cell). Useful for a
        // tile inspector ("who works here?") and for the building-
        // destroyed handler to nudge orphaned agents.
        public IEnumerable<Agent> FindByJob(ChunkCoord coord, int cell)
        {
            for (int i = 0; i < agents.Count; i++)
            {
                Agent a = agents[i];
                if (a == null || !a.Job.HasValue) continue;
                var j = a.Job.Value;
                if (j.Coord.Equals(coord) && j.CellIndex == cell) yield return a;
            }
        }
    }
}
