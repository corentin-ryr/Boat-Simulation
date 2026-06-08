using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using TerrainGrid;

namespace Agents.Actions
{
    // Stand at the assigned Bakery, consume 1 Wheat, produce 1 Bread.
    // Slot-claimed for the duration (one Baker per Bakery in v1).
    //
    // Symbolic contract:
    //   Pre  : AtBakery = true, HasWheat = true, HasBread = false
    //   Eff  : HasWheat = false, HasBread = true
    public sealed class BakeBreadAction : GoapAction
    {
        public const float BakeSeconds = 6f;

        public override string Name => "BakeBread";

        static readonly KeyValuePair<string, object>[] PreCond =
        {
            new KeyValuePair<string, object>(WorldKey.AtBakery, true),
            new KeyValuePair<string, object>(WorldKey.HasWheat, true),
            new KeyValuePair<string, object>(WorldKey.HasBread, false),
        };
        static readonly KeyValuePair<string, object>[] Eff =
        {
            new KeyValuePair<string, object>(WorldKey.HasWheat, false),
            new KeyValuePair<string, object>(WorldKey.HasBread, true),
        };
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => PreCond;
        public override IEnumerable<KeyValuePair<string, object>> Effects => Eff;

        public override float Cost(IGoapAgent agent) => 6f;

        public override bool IsValid(IGoapAgent agent)
        {
            if (!agent.State.Has(TargetKey.AssignedBakery)) return false;
            if (!agent.Inventory.Has(Commodity.Wheat, 1)) return false;
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return false;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedBakery);
            if (!reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e)) return false;
            return e.SlotsClaimed.Contains(agent) || e.HasFreeSlot;
        }

        public override void Start(IGoapAgent agent)
        {
            agent.GetActionData<BakeBreadActionData>().Reset();
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedBakery);
            if (reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
                reg.TryClaimSlot(e, agent);
        }

        public override bool Perform(IGoapAgent agent, float dt)
        {
            BakeBreadActionData data = agent.GetActionData<BakeBreadActionData>();
            World.GameClock clock = World.GameClock.Instance;
            float gameDt = clock != null
                ? dt * (3600f / UnityEngine.Mathf.Max(0.01f, clock.secondsPerGameHour))
                : dt;
            data.ProgressGameSeconds += gameDt;
            return data.ProgressGameSeconds >= BakeSeconds;
        }

        public override void Complete(IGoapAgent agent)
        {
            agent.Inventory.Remove(Commodity.Wheat, 1);
            agent.Inventory.Add(Commodity.Bread, 1);
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedBakery);
            if (reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
                reg.ReleaseSlot(e, agent);
        }

        public override void Stop(IGoapAgent agent)
        {
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.AssignedBakery);
            if (reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
                reg.ReleaseSlot(e, agent);
        }
    }

    public sealed class BakeBreadActionData : ActionData
    {
        public float ProgressGameSeconds;
        public override void Reset() { ProgressGameSeconds = 0f; }
    }
}
