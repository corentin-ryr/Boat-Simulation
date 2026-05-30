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

        public PrimalChunk(ChunkCoord coord, List<Polygon> polygons, VertexCollection verts)
        {
            Coord = coord;
            Polygons = polygons;
            Verts = verts;
            Version = 0;
        }

        // Build a fully independent copy of this chunk's graph (new Vertex/Polygon/VertexCollection
        // objects). Used to publish a chunk to the render-side mirror: the worker keeps mutating
        // the original on later relaxation passes, so the main thread must hold its own copy.
        // Equivalent to ChunkSnapshot.Capture + Restore, but in a single pass with no intermediate
        // primitive arrays. Safe to call off the main thread (pure object construction).
        public PrimalChunk DeepCopy()
        {
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.DeepCopy);
            TerrainProfiler.IncDeepCopies();

            // Pre-size the lookup map and vertex collection so they never rehash. Iterate
            // Verts.Values directly (allocation-free) instead of Verts.ToArray().
            int vCount = Verts.Count;
            Dictionary<Vertex, Vertex> map = new Dictionary<Vertex, Vertex>(vCount);
            VertexCollection vc = new VertexCollection();

            foreach (Vertex v in Verts.Values)
            {
                Vertex nv = new Vertex(v.Position, v.IsEdge);
                map[v] = nv;
                vc.AddVertex(nv);
            }

            List<Polygon> polys = new List<Polygon>(Polygons.Count);
            foreach (Polygon p in Polygons)
            {
                Vertex[] src = p.GetVertices();
                Vertex[] dst = new Vertex[src.Length];
                for (int i = 0; i < src.Length; i++) dst[i] = map[src[i]]; // direct, all keys present
                Polygon np = new Polygon(dst); // ctor re-links vertex -> polygon associations
                np.Terrain = p.Terrain;        // carry the per-cell classification across the copy
                polys.Add(np);
            }

            // Copy the cached dual too. Dual polygons own their own Vertex objects (no
            // shared refs with the primal graph), so we build a fresh independent set —
            // consumer mirrors then never see worker-side mutations to vertex positions.
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
                    dualCopy.Add(new Polygon(dst));
                }
            }

            return new PrimalChunk(Coord, polys, vc)
            {
                Version = Version,
                Dual = dualCopy,
                DualBuiltFromVersion = DualBuiltFromVersion,
                IsFlat = IsFlat,
            };
        }
    }
}
