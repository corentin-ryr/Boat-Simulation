using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using TerrainGrid;

namespace Agents.Actions
{
    // Atomic sale at the agent's nearest Market — 1 unit of `commodity`
    // out of agent.Inventory in exchange for the Market's current price.
    // Instant (Perform returns true on first call) — there is no
    // animation hook in v1.
    //
    // Symbolic effect declares MoneyEnough=true so the planner can chain
    // Harvest → MoveToMarket → Sell to satisfy EarnGoal. Real-world
    // MoneyEnough is projected by MoneySensor from the actual Money
    // field, so if 1 sale doesn't actually cross the threshold the
    // planner just re-plans the same Harvest → Sell loop.
    //
    // Parameterized by Commodity so SellWheat and SellBread are two
    // configured instances in the blueprint.
    public sealed class SellAtMarketAction : GoapAction
    {
        readonly Commodity commodity;
        readonly string actionName;
        readonly KeyValuePair<string, object>[] pre;
        readonly KeyValuePair<string, object>[] eff;

        public SellAtMarketAction(Commodity c)
        {
            commodity = c;
            actionName = "Sell" + Commodities.Name(c);
            string hasKey = c == Commodity.Wheat ? WorldKey.HasWheat : WorldKey.HasBread;
            pre = new[]
            {
                new KeyValuePair<string, object>(WorldKey.AtMarket, true),
                new KeyValuePair<string, object>(hasKey,            true),
            };
            eff = new[]
            {
                new KeyValuePair<string, object>(hasKey,             false),
                new KeyValuePair<string, object>(WorldKey.MoneyEnough, true),
            };
        }

        public override string Name => actionName;
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => pre;
        public override IEnumerable<KeyValuePair<string, object>> Effects => eff;
        public override float Cost(IGoapAgent agent) => 1f;

        public override bool IsValid(IGoapAgent agent)
        {
            if (!agent.Inventory.Has(commodity, 1)) return false;
            return ResolveMarket(agent) != null;
        }

        public override bool Perform(IGoapAgent agent, float dt) => true;

        public override void Complete(IGoapAgent agent)
        {
            MarketStockpile m = ResolveMarket(agent);
            if (m == null) return;
            if (m.TryBuyFrom(commodity, 1, out int paid))
            {
                agent.Inventory.Remove(commodity, 1);
                agent.Money += paid;
            }
        }

        static MarketStockpile ResolveMarket(IGoapAgent agent)
        {
            if (!agent.State.Has(TargetKey.NearestMarket)) return null;
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(TargetKey.NearestMarket);
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return null;
            if (!reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e)) return null;
            return e.Market;
        }
    }
}
