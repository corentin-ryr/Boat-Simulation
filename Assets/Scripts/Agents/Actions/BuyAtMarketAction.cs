using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using TerrainGrid;
using UnityEngine;

namespace Agents.Actions
{
    // Atomic purchase at the agent's nearest Market — 1 unit of
    // `commodity` taken from the Market's Stock in exchange for the
    // current price out of agent.Money. Instant Perform.
    //
    // Symbolic effect just sets HasX=true — we don't model the cash
    // outflow in the planner because Buying isn't pursued by EarnGoal.
    // MoneySensor catches the real Money delta on the next tick.
    public sealed class BuyAtMarketAction : GoapAction
    {
        readonly Commodity commodity;
        readonly string actionName;
        readonly KeyValuePair<string, object>[] pre;
        readonly KeyValuePair<string, object>[] eff;

        public BuyAtMarketAction(Commodity c)
        {
            commodity = c;
            actionName = "Buy" + Commodities.Name(c);
            string hasKey = c == Commodity.Wheat ? WorldKey.HasWheat : WorldKey.HasBread;
            pre = new[]
            {
                new KeyValuePair<string, object>(WorldKey.AtMarket, true),
                new KeyValuePair<string, object>(hasKey,            false),
            };
            eff = new[]
            {
                new KeyValuePair<string, object>(hasKey, true),
            };
        }

        public override string Name => actionName;
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => pre;
        public override IEnumerable<KeyValuePair<string, object>> Effects => eff;
        public override float Cost(IGoapAgent agent) => 1f;

        public override bool IsValid(IGoapAgent agent)
        {
            MarketStockpile m = ResolveMarket(agent);
            if (m == null) return false;
            if (!m.Stock.Has(commodity, 1)) return false;
            int price = Mathf.CeilToInt(m.GetPrice(commodity));
            return agent.Money >= price;
        }

        public override bool Perform(IGoapAgent agent, float dt) => true;

        public override void Complete(IGoapAgent agent)
        {
            MarketStockpile m = ResolveMarket(agent);
            if (m == null) return;
            if (m.TrySellTo(commodity, 1, agent.Money, out int paid))
            {
                agent.Money -= paid;
                agent.Inventory.Add(commodity, 1);
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
