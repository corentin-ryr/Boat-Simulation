using System.Collections.Generic;

namespace Agents.Core
{
    // A goal the agent might pursue. Stateless like GoapAction — a
    // single instance per concrete goal type, evaluated against any
    // agent's current state via Priority(agent).
    //
    // Selection contract (lives in Agent.Tick): pick the highest-
    // priority goal whose TargetState isn't already satisfied. Apply
    // hysteresis (Agent.minGoalSeconds) so a goal once committed to
    // sticks unless its priority drops below a "drop" threshold or
    // the current plan becomes invalid. Without hysteresis, oscillating
    // priorities (Hungry / Tired flipping by 1 each tick) would have
    // the agent thrashing between Eat and Sleep without progressing
    // either.
    public abstract class GoapGoal
    {
        public abstract string Name { get; }

        // The symbolic state the planner must reach for this goal to
        // be considered "satisfied." Same shape as GoapAction.Effects —
        // (key, value) pairs that, when present in the WorldState bag,
        // mean we're done.
        public abstract IEnumerable<KeyValuePair<string, object>> TargetState { get; }

        // Higher = more urgent. Compared cross-goal at goal-selection
        // time. Concrete goals usually project from real-world fields:
        //   SurviveGoal   ∝ max(100 - Hunger, 100 - Energy)
        //   EarnGoal      ∝ max(0, MoneyThreshold - Money)
        //   RestockGoal   ∝ (Bread < threshold) ? fixed : 0
        public abstract float Priority(IGoapAgent agent);

        public override string ToString() => Name;
    }
}
