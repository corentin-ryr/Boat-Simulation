namespace TerrainGrid
{
    // Gameplay classification of a primal cell, set and mutated independently of the
    // generation-time CellTerrain (which is read-only after generation). Drives AI
    // navigation, building placement, and any other "what is this tile?" query.
    public enum TileKind : byte
    {
        Default = 0,
        Road = 1,
        Building = 2,
    }

    // Per-cell gameplay slot. `Payload` is reserved for kind-specific data (e.g. a
    // building id, a road-direction bitmask, a damage counter) — meaning is decided
    // per-kind by gameplay code. Kept small so the parallel array on a chunk stays
    // tight.
    public struct TileProperty
    {
        public TileKind Kind;
        public ushort Payload;
    }

    // Main → worker mutation message. Coord + PolygonIndex identify the target cell
    // (stable across snapshot round-trip because ChunkSnapshot preserves polygon order).
    public struct TileMutation
    {
        public ChunkCoord Coord;
        public int PolygonIndex;
        public TileProperty Value;
    }
}
