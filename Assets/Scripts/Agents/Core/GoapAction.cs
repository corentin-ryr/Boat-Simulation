using System.Collections.Generic;

namespace Agents.Core
{
    // Abstract base for an action the planner can chain into a plan.
    // Stateless contract — no per-agent mutable fields on the action
    // itself; per-agent scratch lives in ActionData buckets on the
    // agent. A single instance per concrete action type is fine for
    // 100 agents in v1.
    //
    // Lifecycle (driven by Agent.Tick when this action is current):
    //   IsValid(agent)   — pre-flight, re-checked every tick. Failure
    //                      invalidates the current plan and forces a
    //                      replan.
    //   Start(agent)     — once when the action becomes current. Reset
    //                      ActionData, claim slots, snap pathfinder.
    //   Perform(agent,dt)
    //                    — every tick. Returns true when the action
    //                      completes (Complete() runs next). Return
    //                      false to keep running.
    //   Complete(agent)  — once on natural completion. Apply real-world
    //                      effects (debit inventory, credit money,
    //                      release slots).
    //   Stop(agent)      — once on premature interruption (plan
    //                      invalidated, agent died, replan). Release
    //                      slots without applying effects.
    //
    // Preconditions / Effects are SYMBOLIC and STATIC — they describe
    // the planner's view of what must be true beforehand and what the
    // action toggles afterwards. The planner mutates a clone of the
    // WorldState by applying Effects; the real-world side runs in
    // Complete().
    public abstract class GoapAction
    {
        public abstract string Name { get; }

        // Symbolic state required for the planner to consider this
        // action. Returned as an IEnumerable so subclasses can either
        // build a List once and cache it, or yield-return.
        public abstract IEnumerable<KeyValuePair<string, object>> Preconditions { get; }

        // Symbolic state mutation the planner applies on expansion.
        public abstract IEnumerable<KeyValuePair<string, object>> Effects { get; }

        // Plan-cost heuristic. Constant 1 by default — concrete actions
        // override to weight long operations (BakeBread) above quick
        // ones (Eat). Cost is summed across the plan; the planner
        // prefers the cheapest total.
        public virtual float Cost(IGoapAgent agent) => 1f;

        // Dynamic guard. Called by the planner to skip impossible
        // actions during expansion AND by Agent.Tick before every
        // Perform() so a plan whose precondition expired (slot taken,
        // inventory exhausted, target destroyed) is dropped.
        //
        // Returning false has different consequences depending on
        // when it fires: during planning ⇒ action skipped that
        // expansion. During execution ⇒ current plan invalidated.
        public virtual bool IsValid(IGoapAgent agent) => true;

        public virtual void Start(IGoapAgent agent) { }

        // The work loop. Return true when done. dt is wall-clock
        // delta (GameClock.timeScale already baked in by Agent).
        public abstract bool Perform(IGoapAgent agent, float dt);

        public virtual void Complete(IGoapAgent agent) { }

        public virtual void Stop(IGoapAgent agent) { }

        public override string ToString() => Name;
    }
}
