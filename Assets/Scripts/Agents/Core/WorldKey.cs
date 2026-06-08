namespace Agents.Core
{
    // Centralised symbolic-state keys for the planner's WorldState bag.
    // Const strings rather than ScriptableObjects — same intent as
    // CrashKonijn's `IWorldKey` types, sized for our code-first design.
    //
    // Adding a key:
    //   1. Add the const here.
    //   2. Add a sensor (or extend an existing one) that writes it.
    //   3. Reference it from any action's Preconditions / Effects.
    //   4. Reference it from any goal's TargetState.
    //
    // Keys ending in `Enough` are sensor-projected booleans for numeric
    // thresholds — they let the planner reason about "is the agent rich
    // enough" without modelling money as a number in the bag.
    public static class WorldKey
    {
        // -------- Clock / time --------
        public const string Hour      = "Hour";       // int 0..23
        public const string IsNight   = "IsNight";    // bool

        // -------- Carrying state (mirrored from Inventory by sensors) --------
        public const string HasWheat  = "HasWheat";   // bool (Wheat >= 1)
        public const string HasBread  = "HasBread";   // bool (Bread >= 1)

        // -------- Needs (mirrored from AgentNeeds) --------
        public const string IsHungry  = "IsHungry";   // bool (Hunger < threshold)
        public const string IsTired   = "IsTired";    // bool (Energy < threshold)
        public const string Rested    = "Rested";     // bool (Energy >= 90)
        public const string Fed       = "Fed";        // bool (Hunger >= 70)

        // -------- Money (mirrored from Agent.Money) --------
        public const string MoneyEnough = "MoneyEnough"; // bool (Money >= goal threshold)

        // -------- Location (where am I right now) --------
        // Stored as (TerrainGrid.ChunkCoord, int) value-tuple. Sensors
        // project Agent's current cell here; MoveToTileAction's effect
        // sets it to the destination cell.
        public const string AtTile = "AtTile";

        // -------- Sub-location flags used by preconditions --------
        // Boolean projections of AtTile against the agent's assigned
        // targets. Sensors set them; actions check them.
        public const string AtField   = "AtField";
        public const string AtBakery  = "AtBakery";
        public const string AtMarket  = "AtMarket";
        public const string AtHome    = "AtHome";
    }
}
