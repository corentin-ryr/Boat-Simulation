using UnityEngine;

namespace TerrainGrid
{
    // Placement rules for paintable TileKinds. Called by TileSelector (RTS
    // click) and TilePainter (gameplay click) before EnqueueTileMutation, so
    // invalid clicks are dropped at the call site rather than being applied
    // and then rolled back.
    //
    // Current rules:
    //   Default / Road / Building / Field / Market — allowed anywhere.
    //   Lighthouse / Dock — only on cells whose primal vertices straddle
    //     sea level (coastal), OR adjacent to an existing Lighthouse/Dock
    //     cell. The adjacency clause lets the player extend a harbour
    //     cluster inland by one ring per click.
    //
    // Adding a new TileKind = extend the switch in IsValid + a new
    // helper method if the rule is novel.
    public static class TilePlacement
    {
        // World Y the ocean surface sits at (boat-sim convention).
        const float SeaLevelY = 0f;
        // Tolerance — a vertex within ±epsilon of sea level counts as
        // neither above nor below. Prevents floating-point noise from
        // pinning a vertex exactly on the ocean plane and failing both
        // halves of the straddle test.
        const float SeaLevelEpsilon = 0.001f;

        public static bool IsValid(TileKind kind, ChunkCoord coord, int cellIdx, ChunkManager mgr)
        {
            switch (kind)
            {
                case TileKind.Lighthouse:
                case TileKind.Dock:
                    return IsCoastalOrAdjacent(coord, cellIdx, mgr, TileKind.Lighthouse, TileKind.Dock);
                default:
                    return true;
            }
        }

        // True if the cell at (coord, cellIdx) is coastal — primal vertices
        // straddle sea level — OR adjacent to a cell of `kindA` / `kindB`,
        // walking the CSR neighbour table and resolving cross-chunk
        // neighbours through the published render mirror.
        static bool IsCoastalOrAdjacent(ChunkCoord coord, int cellIdx, ChunkManager mgr,
                                        TileKind kindA, TileKind kindB)
        {
            if (mgr == null || cellIdx < 0) return false;
            if (!mgr.TryGetPublishedChunk(coord, out PrimalChunk chunk)) return false;

            if (IsCoastalCell(chunk, cellIdx)) return true;
            if (HasNeighborOfKinds(chunk, cellIdx, kindA, kindB, mgr)) return true;
            return false;
        }

        // Walks the CSR edge table of the cell, collecting min/max Y across
        // all edge endpoints. Each primal vertex appears at least once across
        // the cell's edges (twice for interior verts), so the min/max covers
        // every corner.
        static bool IsCoastalCell(PrimalChunk chunk, int cellIdx)
        {
            int[] offsets = chunk.PrimalEdgeOffsets;
            Vector3[] a = chunk.PrimalEdgeVertexA;
            Vector3[] b = chunk.PrimalEdgeVertexB;
            if (offsets == null || a == null || b == null) return false;
            if (cellIdx + 1 >= offsets.Length) return false;

            int start = offsets[cellIdx];
            int end   = offsets[cellIdx + 1];
            if (end <= start) return false;

            bool anyBelow = false;
            bool anyAbove = false;
            for (int j = start; j < end; j++)
            {
                if (j < 0 || j >= a.Length || j >= b.Length) continue;
                if (a[j].y < SeaLevelY - SeaLevelEpsilon) anyBelow = true;
                if (a[j].y > SeaLevelY + SeaLevelEpsilon) anyAbove = true;
                if (b[j].y < SeaLevelY - SeaLevelEpsilon) anyBelow = true;
                if (b[j].y > SeaLevelY + SeaLevelEpsilon) anyAbove = true;
                if (anyBelow && anyAbove) return true;
            }
            return false;
        }

        // True if any primal neighbour of cellIdx (including across chunk
        // seams) carries kindA or kindB. Chunks not currently loaded are
        // skipped — a Lighthouse in an off-screen chunk shouldn't grant
        // placement legality to cells the player can't navigate to anyway.
        static bool HasNeighborOfKinds(PrimalChunk chunk, int cellIdx,
                                       TileKind kindA, TileKind kindB, ChunkManager mgr)
        {
            int[] offsets = chunk.PrimalEdgeOffsets;
            int[] dq      = chunk.PrimalNeighborChunkDQ;
            int[] dr      = chunk.PrimalNeighborChunkDR;
            int[] nbIdx   = chunk.PrimalNeighborPolyIdx;
            if (offsets == null || dq == null || dr == null || nbIdx == null) return false;
            if (cellIdx + 1 >= offsets.Length) return false;

            int start = offsets[cellIdx];
            int end   = offsets[cellIdx + 1];
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
                if (k == kindA || k == kindB) return true;
            }
            return false;
        }
    }
}
