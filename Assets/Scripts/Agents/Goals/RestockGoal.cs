using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;

namespace Agents.Goals
{
    // Keep a stash of Bread in inventory so EatAction has something to
    // consume when hunger creeps up. Plan: MoveToMarket → BuyBread.
    // Bakers naturally restock through their own production chain — for
    // them this goal's chain (BuyBread) competes with Bake. Cost makes
    // BuyBread cheaper, so the planner usually picks the buy unless
    // they're already at the bakery.
    public sealed class RestockGoal : GoapGoal
    {
        public override string Name => "Restock";

        static readonly KeyValuePair<string, object>[] Target =
        {
            new KeyValuePair<string, object>(WorldKey.HasBread, true),
        };
        public override IEnumerable<KeyValuePair<string, object>> TargetState => Target;

        public override float Priority(IGoapAgent agent)
        {
            // Flat 40 when out of bread; quiet otherwise. Survive's
            // urgency curve climbs above 40 once Hunger drops below ~60,
            // so Restock only "wins" when the agent has bread runway.
            return agent.Inventory.Has(Commodity.Bread, 1) ? 0f : 40f;
        }
    }
}
