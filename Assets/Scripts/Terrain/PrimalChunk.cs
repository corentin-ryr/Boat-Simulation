using System.Collections.Generic;

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

            // Primal graph deliberately omitted (Polygons=null, Verts=null). Consumers must
            // null-check if they access them.
            return new PrimalChunk(Coord, null, null)
            {
                Version = Version,
                Dual = dualCopy,
                DualBuiltFromVersion = DualBuiltFromVersion,
                DualComplete = DualComplete,
                IsFlat = IsFlat,
            };
        }
    }
}
