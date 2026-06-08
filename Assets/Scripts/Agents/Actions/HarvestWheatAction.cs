using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using TerrainGrid;

namespace Agents.Actions
{
    // Stand at the assigned Field and harvest 1 Wheat. Slot-claimed for
    // the duration so two farmers can't reap the same furrow. Wheat
    // accumulates on the Field's local Stock (BuildingEntry.Stock); the
    // agent transfers it into their own Inventory in Complete().
    //
    // Symbolic contract:
    //   Pre  : AtField = true, HasWheat = false
    //   Eff  : HasWheat = true
    public sealed class HarvestWheatAction : GoapAction
    {
        // Game-time seconds to complete one harvest. Scales with timeScale
        // via Agent.Update's per-frame Tick.
        public const float HarvestSeconds = 5f;

        // Field's regenerating stockpile cap. Per cell, refilled by an
        // off-screen tick (BuildingRegistry could host this; for v1 we
        // refill greedily on harvest start so the demo never starves
        // out).
        public const int FieldMaxStock = 10;

        public override string Name => "HarvestWheat";

        static readonly KeyValuePair<string, object>[] PreCond =
        {
            new KeyValuePair<string, object>(WorldKey.AtField,  true),
            new KeyValuePair<string, object>(WorldKey.HasWheat, false),
        };
        static readonly KeyValuePair<string, object>[] Eff =
        {
            new KeyValuePair<string, object>(WorldKey.HasWheat, true),
        };
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => PreCond;
        public override IEnumerable<KeyValuePair<string, object>> Effects => Eff;

        public override float Cost(IGoapAgent agent) => 5f;

        public override bool IsValid(IGoapAgent agent)
        {
            if (!agent.State.Has(TargetKey.AssignedField)) return false;
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return false;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedField);
            if (!reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e)) return false;
            // Slot claimable (or already mine) + stock > 0.
            bool slotOk = e.SlotsClaimed.Contains(agent) || e.HasFreeSlot;
            bool stockOk = e.Stock.Has(Commodity.Wheat, 1);
            // v1: lazy regeneration — top up to MaxStock at every IsValid
            // call so the demo isn't gated on a separate refill tick.
            if (!stockOk)
            {
                e.Stock.Set(Commodity.Wheat, FieldMaxStock);
                stockOk = true;
            }
            return slotOk && stockOk;
        }

        public override void Start(IGoapAgent agent)
        {
            agent.GetActionData<HarvestWheatActionData>().Reset();
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedField);
            if (reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
                reg.TryClaimSlot(e, agent);
        }

        public override bool Perform(IGoapAgent agent, float dt)
        {
            HarvestWheatActionData data = agent.GetActionData<HarvestWheatActionData>();
            // dt is wall-time; convert via GameClock so the harvest
            // duration is game-time. Faster timeScale → faster harvest.
            World.GameClock clock = World.GameClock.Instance;
            float gameDt = clock != null
                ? dt * (3600f / UnityEngine.Mathf.Max(0.01f, clock.secondsPerGameHour))
                : dt;
            data.ProgressGameSeconds += gameDt;
            return data.ProgressGameSeconds >= HarvestSeconds;
        }

        public override void Complete(IGoapAgent agent)
        {
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedField);
            if (reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
            {
                e.Stock.Remove(Commodity.Wheat, 1);
                agent.Inventory.Add(Commodity.Wheat, 1);
                reg.ReleaseSlot(e, agent);
            }
        }

        public override void Stop(IGoapAgent agent)
        {
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedField);
            if (reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
                reg.ReleaseSlot(e, agent);
        }
    }

    public sealed class HarvestWheatActionData : ActionData
    {
        public float ProgressGameSeconds;
        public override void Reset() { ProgressGameSeconds = 0f; }
    }
}
