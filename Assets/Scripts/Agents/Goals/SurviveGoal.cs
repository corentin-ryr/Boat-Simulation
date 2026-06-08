using System.Collections.Generic;
using Agents.Core;
using Agents.Runtime;
using UnityEngine;

namespace Agents.Goals
{
    // Keep Hunger and Energy above the IsFed / IsRested thresholds —
    // drives the eat / sleep loop. Priority climbs as either need
    // depletes; a night-time bonus makes Survive easier to pick over
    // Earn once the village's working hours are over.
    public sealed class SurviveGoal : GoapGoal
    {
        public override string Name => "Survive";

        // Both needs satisfied = the goal is "done." MoveToHome / Sleep /
        // (BuyBread + Eat) plan chains close the gap.
        static readonly KeyValuePair<string, object>[] Target =
        {
            new KeyValuePair<string, object>(WorldKey.Fed,    true),
            new KeyValuePair<string, object>(WorldKey.Rested, true),
        };
        public override IEnumerable<KeyValuePair<string, object>> TargetState => Target;

        public override float Priority(IGoapAgent agent)
        {
            Agent concrete = agent as Agent;
            if (concrete == null) return 0f;
            float hungerUrgency = Mathf.Max(0f, 100f - concrete.Needs.Hunger);
            float energyUrgency = Mathf.Max(0f, 100f - concrete.Needs.Energy);
            float nightBonus = agent.State.Get<bool>(WorldKey.IsNight) ? 30f : 0f;
            return Mathf.Max(hungerUrgency, energyUrgency, nightBonus);
        }
    }
}
