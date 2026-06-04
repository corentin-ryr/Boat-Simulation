namespace TerrainGrid
{
    // Thin dispatcher — looks up the TileDefinition for the kind being painted
    // and asks its `placementRule` whether the click is legal. A definition with
    // a null rule (or a kind with no entry in the catalog at all) is treated
    // as "anywhere allowed".
    //
    // Called by TileSelector (RTS click) and TilePainter (gameplay click)
    // before EnqueueTileMutation, so an illegal click is dropped at the call
    // site instead of being applied and then needing back-out logic in the
    // streamer.
    //
    // All actual rule logic (coastal detection, adjacency walks, etc.) lives
    // in the concrete TilePlacementRule subclasses under
    // Assets/Scripts/Terrain/Placement/. To add a new constraint type: create
    // a new TilePlacementRule subclass + asset, drop it onto the relevant
    // TileDefinition's `placementRule` slot.
    public static class TilePlacement
    {
        public static bool IsValid(TileKind kind, ChunkCoord coord, int cellIdx, ChunkManager mgr)
        {
            if (mgr == null) return false;
            TileCatalog catalog = mgr.Catalog;
            if (catalog == null) return true;             // no catalog → no rules → permissive
            TileDefinition def = catalog.Get(kind);
            if (def == null || def.placementRule == null) return true;  // null rule → anywhere
            return def.placementRule.IsValid(coord, cellIdx, mgr);
        }
    }
}
