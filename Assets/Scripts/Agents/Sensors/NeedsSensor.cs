using Agents.Core;
using Agents.Runtime;

namespace Agents.Sensors
{
    // Mirrors AgentNeeds.{Hunger, Energy} into the four threshold-bool
    // WorldKeys the planner reasons about. Casts IGoapAgent down to the
    // concrete Agent because Needs is a Runtime concept — Core stays
    // independent of Runtime via this layering.
    public sealed class NeedsSensor : GoapSensor
    {
        public override string Name => "Needs";
        public override SenseMode Mode => SenseMode.Always;

        public override void Sense(IGoapAgent agent)
        {
            Agent concrete = agent as Agent;
            if (concrete == null) return;
            AgentNeeds n = concrete.Needs;
            agent.State.Set(WorldKey.IsHungry, n.IsHungry());
            agent.State.Set(WorldKey.IsTired,  n.IsTired());
            agent.State.Set(WorldKey.Fed,      n.IsFed());
            agent.State.Set(WorldKey.Rested,   n.IsRested());
        }
    }
}
