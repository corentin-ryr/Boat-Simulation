namespace Agents.Core
{
    // Spatial-target keys carried by the agent. Stored in WorldState as
    // (TerrainGrid.ChunkCoord, int) value tuples — same shape as a
    // BuildingEntry handle so actions can hand them straight to the
    // pathfinder. Distinct from WorldKey so sensors and actions can tell
    // "what is true?" (WorldKey) apart from "where is it?" (TargetKey).
    //
    // Adding a target:
    //   1. Add the const here.
    //   2. Add a sensor that resolves a BuildingEntry via
    //      BuildingRegistry.FindNearestByRole and stores its
    //      (Coord, CellIndex) tuple.
    //   3. Reference it from MoveToTileAction's destination resolution.
    public static class TargetKey
    {
        // The Farmer's "I work this Field" assignment. Stored on spawn by
        // NearestEmploymentSensor and refreshed when the Field is destroyed.
        public const string AssignedField  = "AssignedField";

        // The Baker's workshop assignment.
        public const string AssignedBakery = "AssignedBakery";

        // Nearest Market (re-resolved as the agent moves, since "nearest"
        // depends on the agent's current position).
        public const string NearestMarket  = "NearestMarket";

        // The agent's residence.
        public const string Home           = "Home";
    }
}
