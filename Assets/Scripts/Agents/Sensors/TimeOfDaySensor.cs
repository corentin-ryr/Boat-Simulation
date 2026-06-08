using Agents.Core;
using World;

namespace Agents.Sensors
{
    // Projects GameClock.GameHour + day/night flag into the agent's
    // WorldState. Interval mode — checking every game-minute is plenty
    // fast for goal evaluation; planners don't care about sub-minute
    // resolution.
    public sealed class TimeOfDaySensor : GoapSensor
    {
        public override string Name => "TimeOfDay";
        public override SenseMode Mode => SenseMode.Interval;

        // Game-seconds between updates. A game-minute = 60 game-seconds.
        public override float IntervalSeconds => 60f;

        public override void Sense(IGoapAgent agent)
        {
            GameClock clock = GameClock.Instance;
            if (clock == null) return;
            agent.State.Set(WorldKey.Hour, clock.GameHour);
            agent.State.Set(WorldKey.IsNight, clock.IsNighttime);
        }
    }
}
