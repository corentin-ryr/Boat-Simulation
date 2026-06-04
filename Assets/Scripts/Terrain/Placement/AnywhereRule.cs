using UnityEngine;

namespace TerrainGrid.Placement
{
    // Permissive rule — accepts every cell. Mostly useful as the explicit
    // default when you want to make the "no restrictions" call site read
    // symmetrically with the constrained kinds. A null `placementRule` on
    // TileDefinition is treated identically (anywhere is the default).
    [CreateAssetMenu(menuName = "Terrain/Placement Rule/Anywhere", fileName = "AnywhereRule")]
    public class AnywhereRule : TilePlacementRule
    {
        public override bool IsValid(ChunkCoord coord, int cellIdx, ChunkManager mgr) => true;
    }
}
