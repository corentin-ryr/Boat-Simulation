namespace TerrainGrid
{
    // Gameplay classification of a primal cell, set and mutated independently of the
    // generation-time CellTerrain (which is read-only after generation). Drives AI
    // navigation, building placement, and any other "what is this tile?" query.
    //
    // Carving (terrain wedge skipped → overlay paints the cell):
    //   Default — not carved (terrain shows through).
    //   Road, Building, Field — carved; ChunkSurface.BuildMesh skips the wedge and the
    //     kind's overlay paints the cell pixel-for-pixel.
    //   Market, Lighthouse, Dock — not carved; the structure mesh sits on top of the
    //     normal terrain wedge so the player sees ground around / under the structure.
    //
    // Placement (TilePlacement.IsValid):
    //   Default, Road, Building, Field, Market — paintable anywhere.
    //   Lighthouse, Dock — only on coastal cells (primal vertices straddle sea level)
    //     OR adjacent to an existing Lighthouse/Dock cell. The adjacency clause lets
    //     the player extend a small harbour cluster inland by one ring per click.
    public enum TileKind : byte
    {
        Default    = 0,
        Road       = 1,
        Building   = 2,
        Field      = 3,
        Market     = 4,
        Lighthouse = 5,
        Dock       = 6,
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
