using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using Agents.Runtime;

namespace Agents.Actions
{
    // Consume 1 Bread. Refills the agent's Hunger by Needs.HungerPerBread.
    // No location requirement — eat wherever you are.
    //
    // Symbolic contract:
    //   Pre  : HasBread = true, Fed = false
    //   Eff  : HasBread = false, Fed = true
    public sealed class EatAction : GoapAction
    {
        public override string Name => "Eat";

        static readonly KeyValuePair<string, object>[] PreCond =
        {
            new KeyValuePair<string, object>(WorldKey.HasBread, true),
            new KeyValuePair<string, object>(WorldKey.Fed,      false),
        };
        static readonly KeyValuePair<string, object>[] Eff =
        {
            new KeyValuePair<string, object>(WorldKey.HasBread, false),
            new KeyValuePair<string, object>(WorldKey.Fed,      true),
        };
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => PreCond;
        public override IEnumerable<KeyValuePair<string, object>> Effects => Eff;
        public override float Cost(IGoapAgent agent) => 1f;

        public override bool IsValid(IGoapAgent agent) => agent.Inventory.Has(Commodity.Bread, 1);

        public override bool Perform(IGoapAgent agent, float dt) => true;

        public override void Complete(IGoapAgent agent)
        {
            agent.Inventory.Remove(Commodity.Bread, 1);
            Agent concrete = agent as Agent;
            if (concrete != null) concrete.Needs.EatBread();
        }
    }
}
