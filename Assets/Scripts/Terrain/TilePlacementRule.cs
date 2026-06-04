using UnityEngine;

namespace TerrainGrid
{
    // A placement rule decides whether the user is allowed to paint a tile of a
    // given kind onto a specific cell. ScriptableObject so the same C# class
    // can be reused with different parameters — e.g. CoastalOrAdjacentRule
    // configured with {Lighthouse, Dock} for harbour structures vs.
    // {Building, Road} for a future "Plaza" rule.
    //
    // Rules ONLY read from the published render mirror (via ChunkManager
    // accessors). They must not enqueue mutations, mutate state, or block.
    // Called per click — so cheap, even at 60+ FPS during a drag-paint.
    public abstract class TilePlacementRule : ScriptableObject
    {
        public abstract bool IsValid(ChunkCoord coord, int cellIdx, ChunkManager mgr);
    }
}
