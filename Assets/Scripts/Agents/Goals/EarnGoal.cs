using System.Collections.Generic;
using Agents.Core;
using UnityEngine;

namespace Agents.Goals
{
    // Top up the agent's Money to MoneyEnoughThreshold so they can buy
    // bread later. Plan chains depend on what job the agent has:
    //   Farmer  : MoveToField → Harvest → MoveToMarket → SellWheat
    //   Baker   : MoveToMarket → BuyWheat → MoveToBakery → Bake →
    //             MoveToMarket → SellBread
    //   Citizen : no chain reaches MoneyEnough → priority is 0 once Job
    //             is unset.
    public sealed class EarnGoal : GoapGoal
    {
        public override string Name => "Earn";

        static readonly KeyValuePair<string, object>[] Target =
        {
            new KeyValuePair<string, object>(WorldKey.MoneyEnough, true),
        };
        public override IEnumerable<KeyValuePair<string, object>> TargetState => Target;

        public override float Priority(IGoapAgent agent)
        {
            // No job → can't earn. Sell actions still exist but the
            // planner won't find a chain that reaches MoneyEnough
            // without going through a production action.
            if (agent.Job == null) return 0f;
            float urgency = Mathf.Max(0f, 60f - agent.Money);
            return urgency;
        }
    }
}
