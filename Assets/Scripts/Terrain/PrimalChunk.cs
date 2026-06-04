using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid
{
    // One chunk's primal (data-layer) grid. Pure data — no rendering, no Unity objects.
    public class PrimalChunk
    {
        public readonly ChunkCoord Coord;
        public List<Polygon> Polygons;
        public VertexCollection Verts;

        // Incremented when this chunk's visible state changes — either border-relaxation
        // moved its vertices (primal change) OR a neighbour cascade invalidated this
        // chunk's cached dual (dual change with no primal change). The streamer's publish
        // loop uses Version to detect what to re-send to each consumer.
        public int Version;

        // Cached dual: this chunk's owned cells (under deterministic "lowest ChunkCoord
        // wins" ownership). Built by the streamer worker after relaxation and before
        // publish, so consumers receive it ready-made — the main thread doesn't compute
        // polygon math. Stale when DualBuiltFromVersion != Version; invalidated either by
        // the chunk's own primal moving or by a neighbour's primal moving (the cascade is
        // performed in TerrainModel — see InvalidateNeighborDuals).
        public List<Polygon> Dual;
        public int DualBuiltFromVersion = -1;

        // True iff every one of this chunk's 6 neighbour coords was loaded when this Dual
        // was built. When true, the dual already accounts for any of those neighbours that
        // get added later (cache revive, store restore, fresh generate) — ownership and
        // seam cells were resolved with full information — so the add-cascade can skip
        // invalidating this chunk. Reset to false on InvalidateDual; recomputed in BuildDual.
        public bool DualComplete;

        // True iff every primal vertex in this chunk sits at Y=0 — i.e. every noise sample
        // came back at or below the elevation threshold. Set once at generation. When true the
        // streamer skips dual computation entirely and the surface mounts a shared flat-ocean
        // tile instead of building a per-cell mesh — the dominant optimization for open ocean.
        public bool IsFlat;

        // Gameplay slot per primal cell, parallel to Polygons by index. Lazily allocated:
        // null means "every cell is TileProperty.Default" so fresh, never-painted chunks
        // pay zero memory. Worker-thread mutation only (see TerrainModel.ApplyTileMutation
        // / TerrainStreamer.EnqueueTileMutation); main thread reads through the streamer's
        // lock-protected getter.
        public TileProperty[] TileProperties;

        // Cached primal cell centroids, parallel to Polygons by index. Computed once at
        // generation (and restored from snapshot) so the main-thread TilePainter can do a
        // nearest-centroid lookup against a hit world position without having to walk the
        // primal polygon graph (which DeepCopy drops on the published mirror). Shipped
        // through DeepCopy so consumers see it.
        public Vector3[] PrimalCellCentroids;

        // Per dual-mesh-vertex primal cell index, flattened in the same order
        // ChunkSurface.BuildMesh walks (each dual polygon contributes its corners in
        // post-sort order). -1 marks a corner whose owning primal face lives in a neighbour
        // chunk (border completion). Lets the render side palette-map vertex colours from
        // TileProperties without having to reconstruct the primal→dual centroid mapping.
        // Null on flat chunks (they render the shared flat tile, no per-vertex tint).
        public int[] DualVertexPrimalIndex;

        // --- Primal neighbour table (for road / building mesh construction) ---
        //
        // CSR layout describing, for each primal cell, every edge's two endpoint positions
        // and which other primal cell (own chunk OR a neighbour chunk) lies across it.
        // Built by TerrainModel.BuildPrimalNeighborTable after BuildDual / restore / add so
        // cross-chunk neighbours are resolved using the same context the dual builder uses
        // (LatticeKey lookup against currently loaded neighbour chunks).
        //
        // Shipped through DeepCopy onto the published mirror so the main thread can build
        // road geometry (centroid → edge midpoint arms, connection bitmask from neighbouring
        // chunks' TileProperties) without ever touching the worker-side primal graph (which
        // DeepCopy drops, see comment near the bottom of this file).
        //
        // Index conventions per edge slot j ∈ [PrimalEdgeOffsets[i], PrimalEdgeOffsets[i+1]):
        //   PrimalEdgeVertexA[j], PrimalEdgeVertexB[j]   — the two world-space endpoints
        //     of the edge (mid = (A+B)/2, length = |B-A|, perp = (B-A)/length). Both
        //     chunks sharing a seam edge see the same two positions, which is what makes
        //     road arms meet exactly at the midpoint with matching widths.
        //   (PrimalNeighborChunkDQ[j], PrimalNeighborChunkDR[j]) — axial offset of the
        //     neighbouring chunk relative to this chunk's coord. (0,0) = same chunk.
        //   PrimalNeighborPolyIdx[j] — index of the neighbour cell inside that chunk's
        //     Polygons / TileProperties / PrimalCellCentroids arrays. -1 means the edge
        //     is on the boundary of the loaded set (no neighbour chunk loaded at build
        //     time, or chunk-seam edge where the neighbour couldn't be resolved). The
        //     dual-complete cascade in BuildDual refreshes this when a neighbour appears.
        //
        // Null on flat chunks (they don't paint roads — they mount the shared flat tile).
        public int[]     PrimalEdgeOffsets;
        public Vector3[] PrimalEdgeVertexA;
        public Vector3[] PrimalEdgeVertexB;
        public int[]     PrimalNeighborChunkDQ;
        public int[]     PrimalNeighborChunkDR;
        public int[]     PrimalNeighborPolyIdx;

        public PrimalChunk(ChunkCoord coord, List<Polygon> polygons, VertexCollection verts)
        {
            Coord = coord;
            Polygons = polygons;
            Verts = verts;
            Version = 0;
        }

        // Build a fully independent *render-mirror* copy of this chunk: independent dual
        // polygon graph (so the worker can keep mutating positions on its own copy without
        // racing the renderer), but NO primal graph — the consumer never reads primal
        // Polygons / Verts on the published copy. The only consumer-side touch of the primal
        // was ChunkManager's `showPrimalGizmos` debug branch; it now gracefully no-ops on
        // mirrored chunks (worker-side primal is off the main thread anyway).
        //
        // Why this matters: the primal copy used to dominate DeepCopy (one Dictionary, one
        // VertexCollection with internal dict, ~100 Vertex objects each carrying a HashSet,
        // ~50 Polygons + Vertex[] arrays, ~150 vertex→polygon HashSet.Adds). All of that is
        // gone. The dual block also uses Polygon.CreateUnlinked to skip wiring back-pointers
        // the consumer never walks, and lazy Vertex.polygons means the unlinked dual vertices
        // never allocate a HashSet at all.
        //
        // Safe to call off the main thread (pure object construction).
        public PrimalChunk DeepCopy()
        {
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.DeepCopy);
            TerrainProfiler.IncDeepCopies();

            // Dual: dual polygons own their own Vertex objects (no shared refs with the
            // primal graph), so we build a fresh independent set — consumer mirrors then
            // never see worker-side mutations. CreateUnlinked skips the vertex→polygon
            // back-pointer wiring the consumer doesn't need.
            List<Polygon> dualCopy = null;
            if (Dual != null)
            {
                int dCount = Dual.Count;
                dualCopy = new List<Polygon>(dCount);
                for (int pi = 0; pi < dCount; pi++)
                {
                    Polygon p = Dual[pi];
                    Vertex[] src = p.GetVertices();
                    int n = src.Length;
                    Vertex[] dst = new Vertex[n];
                    for (int i = 0; i < n; i++)
                        dst[i] = new Vertex(src[i].Position, src[i].IsEdge);
                    dualCopy.Add(Polygon.CreateUnlinked(dst));
                }
            }

            // Render-side parallel arrays. All three are blittable, so each is a single
            // Array.Copy with no per-element allocation. Centroids + dual index are needed
            // by TilePainter (cursor → tile lookup) and BuildMesh (vertex-colour palette);
            // TileProperties are needed by BuildMesh to choose the colour. Skipped on flat
            // chunks because the shared flat tile path doesn't read them.
            Vector3[] centroidsCopy = null;
            if (PrimalCellCentroids != null)
            {
                centroidsCopy = new Vector3[PrimalCellCentroids.Length];
                System.Array.Copy(PrimalCellCentroids, centroidsCopy, PrimalCellCentroids.Length);
            }
            int[] dvpiCopy = null;
            if (DualVertexPrimalIndex != null)
            {
                dvpiCopy = new int[DualVertexPrimalIndex.Length];
                System.Array.Copy(DualVertexPrimalIndex, dvpiCopy, DualVertexPrimalIndex.Length);
            }
            TileProperty[] tpCopy = null;
            if (TileProperties != null)
            {
                tpCopy = new TileProperty[TileProperties.Length];
                System.Array.Copy(TileProperties, tpCopy, TileProperties.Length);
            }

            // Primal neighbour table — six blittable arrays (CSR offsets + two endpoint
            // streams + three neighbour-reference streams). The main thread reads these to
            // build road geometry; they share the same lifecycle as the other render arrays
            // (refreshed every time the worker rebuilds the dual / restores / adds the chunk).
            int[]     edgeOffCopy = null;
            Vector3[] edgeACopy = null;
            Vector3[] edgeBCopy = null;
            int[]     neighDQCopy = null;
            int[]     neighDRCopy = null;
            int[]     neighIdxCopy = null;
            if (PrimalEdgeOffsets != null)
            {
                edgeOffCopy = new int[PrimalEdgeOffsets.Length];
                System.Array.Copy(PrimalEdgeOffsets, edgeOffCopy, PrimalEdgeOffsets.Length);
            }
            if (PrimalEdgeVertexA != null)
            {
                edgeACopy = new Vector3[PrimalEdgeVertexA.Length];
                System.Array.Copy(PrimalEdgeVertexA, edgeACopy, PrimalEdgeVertexA.Length);
            }
            if (PrimalEdgeVertexB != null)
            {
                edgeBCopy = new Vector3[PrimalEdgeVertexB.Length];
                System.Array.Copy(PrimalEdgeVertexB, edgeBCopy, PrimalEdgeVertexB.Length);
            }
            if (PrimalNeighborChunkDQ != null)
            {
                neighDQCopy = new int[PrimalNeighborChunkDQ.Length];
                System.Array.Copy(PrimalNeighborChunkDQ, neighDQCopy, PrimalNeighborChunkDQ.Length);
            }
            if (PrimalNeighborChunkDR != null)
            {
                neighDRCopy = new int[PrimalNeighborChunkDR.Length];
                System.Array.Copy(PrimalNeighborChunkDR, neighDRCopy, PrimalNeighborChunkDR.Length);
            }
            if (PrimalNeighborPolyIdx != null)
            {
                neighIdxCopy = new int[PrimalNeighborPolyIdx.Length];
                System.Array.Copy(PrimalNeighborPolyIdx, neighIdxCopy, PrimalNeighborPolyIdx.Length);
            }

            // Primal graph deliberately omitted (Polygons=null, Verts=null). Consumers must
            // null-check if they access them.
            return new PrimalChunk(Coord, null, null)
            {
                Version = Version,
                Dual = dualCopy,
                DualBuiltFromVersion = DualBuiltFromVersion,
                DualComplete = DualComplete,
                IsFlat = IsFlat,
                PrimalCellCentroids = centroidsCopy,
                DualVertexPrimalIndex = dvpiCopy,
                TileProperties = tpCopy,
                PrimalEdgeOffsets = edgeOffCopy,
                PrimalEdgeVertexA = edgeACopy,
                PrimalEdgeVertexB = edgeBCopy,
                PrimalNeighborChunkDQ = neighDQCopy,
                PrimalNeighborChunkDR = neighDRCopy,
                PrimalNeighborPolyIdx = neighIdxCopy,
            };
        }
    }
}
