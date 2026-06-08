namespace Agents.Core
{
    // Sensor: bridge from the agent's real-world state (Position,
    // Inventory, Money, AtTile) to the symbolic WorldState bag the
    // planner reads. Stateless like everything else in Core — one
    // instance per concrete sensor, registered globally on the agent
    // runtime, tickled per Mode against every agent.
    //
    // Run modes (Agent.TickSensors enforces these):
    //   Always   — every plan tick. Cheapest sensors only (mirror an
    //              int field into a bool flag), so cost stays bounded
    //              even at 100 agents × 4 Hz.
    //   Once     — first plan tick after spawn / after the goal/job
    //              changes. NearestEmploymentSensor uses this — work
    //              assignment is sticky once chosen.
    //   Interval — every `IntervalSeconds` of game time. For sensors
    //              whose source data changes slowly (Hour-of-day,
    //              Market prices).
    //
    // Event-driven variants override OnEvent — Phase F doesn't ship
    // them yet but the API is here so we don't carve a poll-only path
    // through every sensor.
    public abstract class GoapSensor
    {
        public enum SenseMode : byte { Always, Once, Interval }

        public abstract string Name { get; }
        public abstract SenseMode Mode { get; }

        // Only consulted when Mode == Interval. Game-seconds between
        // re-evaluations.
        public virtual float IntervalSeconds => 1f;

        // Mirror real-world state into agent.State. The sensor reads
        // whatever it needs from the agent (Position, Inventory,
        // BuildingRegistry, GameClock, …) and writes the resulting
        // keys via agent.State.Set.
        public abstract void Sense(IGoapAgent agent);

        public override string ToString() => Name;
    }
}
