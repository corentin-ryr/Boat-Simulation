using System.Collections.Generic;
using Agents.Core;
using Agents.Pathfinding;
using TerrainGrid;
using UnityEngine;

namespace Agents.Actions
{
    // Walk the agent from its current tile to the cell stored in a
    // TargetKey. Parameterized so one class covers MoveToField,
    // MoveToBakery, MoveToMarket, MoveToHome — the village blueprint
    // instantiates four configured copies.
    //
    // Symbolic contract:
    //   Pre  : atKey = false        (we're not already there)
    //   Eff  : atKey = true, all other At-flags = false
    //
    // Real-world contract: tween transform.position along a TilePath
    // resolved at Start(). One ActionData bucket shared across the four
    // instances — only one MoveTo runs on an agent at a time, so the
    // bucket is safe to share. Path lookups go through TilePathCache so
    // adjacent agents heading to the same field reuse the planned route.
    public sealed class MoveToTileAction : GoapAction
    {
        readonly string actionName;
        readonly string atKey;
        readonly string targetKey;
        readonly KeyValuePair<string, object>[] preconditions;
        readonly KeyValuePair<string, object>[] effects;
        readonly ChunkManager mgr;
        readonly TilePathCache cache;

        // Wall-clock metres per second. Game-time decoupled: at faster
        // timeScale, agents still move at the same visible speed —
        // animation looks normal regardless of how much in-game work is
        // crammed into a second.
        public float WalkSpeed = 4f;

        public override string Name => actionName;
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => preconditions;
        public override IEnumerable<KeyValuePair<string, object>> Effects => effects;

        public MoveToTileAction(string name, string atKey, string targetKey,
                                string[] otherAtKeysToClear,
                                ChunkManager mgr, TilePathCache cache)
        {
            this.actionName = name;
            this.atKey = atKey;
            this.targetKey = targetKey;
            this.mgr = mgr;
            this.cache = cache;

            preconditions = new[]
            {
                new KeyValuePair<string, object>(atKey, false),
            };

            List<KeyValuePair<string, object>> eff = new List<KeyValuePair<string, object>>
            {
                new KeyValuePair<string, object>(atKey, true),
            };
            if (otherAtKeysToClear != null)
                foreach (string ok in otherAtKeysToClear)
                    eff.Add(new KeyValuePair<string, object>(ok, false));
            effects = eff.ToArray();
        }

        public override bool IsValid(IGoapAgent agent)
        {
            return agent.State.Has(targetKey);
        }

        public override void Start(IGoapAgent agent)
        {
            MoveToTileActionData data = agent.GetActionData<MoveToTileActionData>();
            data.Reset();
            if (!agent.State.Has(targetKey)) { data.Path = null; return; }
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(targetKey);
            var from = ResolveCurrentCell(agent);
            if (from.CellIndex < 0) { data.Path = null; return; }
            data.Path = cache.GetOrFind(from.Coord, from.CellIndex, tup.Coord, tup.CellIndex);
            data.StepIndex = 0;
        }

        public override bool Perform(IGoapAgent agent, float dt)
        {
            MoveToTileActionData data = agent.GetActionData<MoveToTileActionData>();
            // Path resolution failure (target chunk unloaded mid-execution,
            // or unreachable). Complete-as-done so the plan teardown runs;
            // the next plan tick's sensor will see AtX is still false and
            // re-plan with fresh assumptions.
            if (data.Path == null || data.Path.Length == 0) return true;
            if (data.StepIndex >= data.Path.Length) return true;

            Vector3 stepPos = data.Path[data.StepIndex].WorldPos;
            Vector3 cur = agent.Position;
            Vector3 dir = stepPos - cur;
            dir.y = 0f;
            float distSq = dir.sqrMagnitude;

            // Step arrival threshold — once close enough to the next
            // centroid, hop to the following one. Squared compare to
            // skip the sqrt.
            const float StepArriveR = 0.35f;
            if (distSq < StepArriveR * StepArriveR)
            {
                data.StepIndex++;
                return data.StepIndex >= data.Path.Length;
            }

            float dist = Mathf.Sqrt(distSq);
            float advance = WalkSpeed * dt;
            if (advance > dist) advance = dist;
            Vector3 move = dir / dist * advance;
            Vector3 newPos = cur + move;
            newPos.y = stepPos.y;     // snap to terrain Y along the path
            agent.Position = newPos;
            return false;
        }

        public override void Complete(IGoapAgent agent)
        {
            if (!agent.State.Has(targetKey)) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(targetKey);
            agent.AtTile = (tup.Coord, tup.CellIndex);
            // Symbolic At-flags get set by sensors on the next tick.
            // The planner's clone of WorldState already applied them
            // during the search; the live state catches up naturally.
        }

        // Nearest-centroid lookup in the chunk under the agent's XZ —
        // same pattern as TilePicker. Used at Start to feed the
        // pathfinder a valid "from" cell.
        (ChunkCoord Coord, int CellIndex) ResolveCurrentCell(IGoapAgent agent)
        {
            Vector3 pos = agent.Position;
            ChunkCoord coord = ChunkCoord.FromWorldPos(pos, mgr.HexRadius, mgr.ChunkGridSize);
            if (!mgr.TryGetPublishedChunk(coord, out PrimalChunk chunk)) return (coord, -1);
            Vector3[] cs = chunk.PrimalCellCentroids;
            if (cs == null || cs.Length == 0) return (coord, -1);
            int best = -1;
            float bestSq = float.PositiveInfinity;
            for (int i = 0; i < cs.Length; i++)
            {
                float dx = cs[i].x - pos.x;
                float dz = cs[i].z - pos.z;
                float sq = dx * dx + dz * dz;
                if (sq < bestSq) { bestSq = sq; best = i; }
            }
            return (coord, best);
        }
    }

    // Shared bucket across all MoveTo instances. Only one runs at a
    // time, so this is safe — each Start() resets it.
    public sealed class MoveToTileActionData : ActionData
    {
        public TilePath Path;
        public int StepIndex;
        public override void Reset() { Path = null; StepIndex = 0; }
    }
}
