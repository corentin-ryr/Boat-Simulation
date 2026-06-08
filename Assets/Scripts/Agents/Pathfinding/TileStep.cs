using TerrainGrid;
using UnityEngine;

namespace Agents.Pathfinding
{
    // One node on a TilePath: the primal-cell identifier and a cached
    // world-space centroid. The centroid is captured at FindPath time —
    // chunk reloads later may invalidate the position field, but the
    // (Coord, CellIndex) pair stays addressable. Agents revalidate via
    // BuildingRegistry / ChunkManager whenever they actually consume a
    // step.
    public readonly struct TileStep
    {
        public readonly ChunkCoord Coord;
        public readonly int CellIndex;
        public readonly Vector3 WorldPos;

        public TileStep(ChunkCoord coord, int cellIndex, Vector3 worldPos)
        {
            Coord = coord;
            CellIndex = cellIndex;
            WorldPos = worldPos;
        }

        public override string ToString() => $"{Coord}:{CellIndex}";
    }
}
