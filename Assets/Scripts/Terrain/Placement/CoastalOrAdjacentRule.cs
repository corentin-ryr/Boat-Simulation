using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Placement
{
    // Cell qualifies if either:
    //   1. It's coastal — at least one primal vertex strictly below sea level AND
    //      at least one strictly above, so the cell straddles the shoreline.
    //   2. A primal neighbour (within own chunk or across a chunk seam) carries
    //      one of `adjacentKinds`. This is the "extend the harbour inland by one
    //      ring per click" affordance.
    //
    // Reusable across kinds: assign one CoastalOrAdjacentRule asset to both
    // Lighthouse and Dock with adjacentKinds = {Lighthouse, Dock}, and the rule
    // behaves identically for both. A future Pier kind with adjacentKinds =
    // {Dock} could reuse the C# class via a separate asset.
    [CreateAssetMenu(menuName = "Terrain/Placement Rule/Coastal Or Adjacent", fileName = "CoastalOrAdjacentRule")]
    public class CoastalOrAdjacentRule : TilePlacementRule
    {
        [Tooltip("World Y the world ocean sits at. Cells straddling this plane (≥1 vertex " +
                 "above, ≥1 below) count as coastal.")]
        public float seaLevelY = 0f;

        [Tooltip("Tolerance around seaLevelY. A vertex within ±epsilon counts as neither " +
                 "above nor below — protects against noise pinning a vertex exactly on the " +
                 "ocean plane.")]
        public float seaLevelEpsilon = 0.001f;

        [Tooltip("Cells whose primal neighbours have any of these kinds also qualify. Empty " +
                 "list = coastal cells only.")]
        public List<TileKind> adjacentKinds = new List<TileKind>();

        public override bool IsValid(ChunkCoord coord, int cellIdx, ChunkManager mgr)
        {
            if (mgr == null || cellIdx < 0) return false;
            if (!mgr.TryGetPublishedChunk(coord, out PrimalChunk chunk)) return false;

            if (IsCoastalCell(chunk, cellIdx)) return true;
            if (adjacentKinds.Count > 0 && HasAdjacentOfKinds(chunk, cellIdx, mgr)) return true;
            return false;
        }

        // Walks the CSR edge table of the cell, collecting min/max Y across all
        // edge endpoints. Both endpoints of every edge contribute, so each vertex
        // is sampled at least once (twice for interior verts shared by two edges).
        bool IsCoastalCell(PrimalChunk chunk, int cellIdx)
        {
            int[] offsets = chunk.PrimalEdgeOffsets;
            Vector3[] a = chunk.PrimalEdgeVertexA;
            Vector3[] b = chunk.PrimalEdgeVertexB;
            if (offsets == null || a == null || b == null) return false;
            if (cellIdx + 1 >= offsets.Length) return false;

            int start = offsets[cellIdx];
            int end = offsets[cellIdx + 1];
            if (end <= start) return false;

            bool anyBelow = false, anyAbove = false;
            for (int j = start; j < end; j++)
            {
                if (j < 0 || j >= a.Length || j >= b.Length) continue;
                if (a[j].y < seaLevelY - seaLevelEpsilon) anyBelow = true;
                if (a[j].y > seaLevelY + seaLevelEpsilon) anyAbove = true;
                if (b[j].y < seaLevelY - seaLevelEpsilon) anyBelow = true;
                if (b[j].y > seaLevelY + seaLevelEpsilon) anyAbove = true;
                if (anyBelow && anyAbove) return true;
            }
            return false;
        }

        // Resolves cross-chunk neighbours through the published render mirror.
        // Chunks not currently loaded are skipped — a Lighthouse in an off-screen
        // chunk shouldn't grant placement legality to cells the player can't see
        // anyway.
        bool HasAdjacentOfKinds(PrimalChunk chunk, int cellIdx, ChunkManager mgr)
        {
            int[] offsets = chunk.PrimalEdgeOffsets;
            int[] dq = chunk.PrimalNeighborChunkDQ;
            int[] dr = chunk.PrimalNeighborChunkDR;
            int[] nbIdx = chunk.PrimalNeighborPolyIdx;
            if (offsets == null || dq == null || dr == null || nbIdx == null) return false;
            if (cellIdx + 1 >= offsets.Length) return false;

            int start = offsets[cellIdx];
            int end = offsets[cellIdx + 1];
            for (int j = start; j < end; j++)
            {
                if (j < 0 || j >= dq.Length) continue;
                int nIdx = nbIdx[j];
                if (nIdx < 0) continue;

                ChunkCoord nCoord = new ChunkCoord(chunk.Coord.q + dq[j], chunk.Coord.r + dr[j]);
                PrimalChunk nChunk;
                if (dq[j] == 0 && dr[j] == 0) nChunk = chunk;
                else if (!mgr.TryGetPublishedChunk(nCoord, out nChunk)) continue;

                TileProperty[] nTp = nChunk.TileProperties;
                if (nTp == null || nIdx >= nTp.Length) continue;
                TileKind k = nTp[nIdx].Kind;
                for (int a = 0; a < adjacentKinds.Count; a++)
                    if (adjacentKinds[a] == k) return true;
            }
            return false;
        }
    }
}
