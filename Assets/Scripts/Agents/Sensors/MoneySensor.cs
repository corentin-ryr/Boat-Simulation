using Agents.Core;

namespace Agents.Sensors
{
    // Projects Agent.Money into the symbolic MoneyEnough bool. The
    // threshold matches EarnGoal's TargetState — once Money reaches
    // it, EarnGoal becomes "satisfied" and other goals can take over.
    public sealed class MoneySensor : GoapSensor
    {
        public const int MoneyEnoughThreshold = 60;

        public override string Name => "Money";
        public override SenseMode Mode => SenseMode.Always;

        public override void Sense(IGoapAgent agent)
        {
            agent.State.Set(WorldKey.MoneyEnough, agent.Money >= MoneyEnoughThreshold);
        }
    }
}
