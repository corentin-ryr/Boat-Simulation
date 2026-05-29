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
            Dictionary<Vertex, Vertex> map = new Dictionary<Vertex, Vertex>();
            VertexCollection vc = new VertexCollection();

            Vertex Map(Vertex v)
            {
                if (!map.TryGetValue(v, out Vertex nv))
                {
                    nv = new Vertex(v.Position, v.IsEdge);
                    map[v] = nv;
                    vc.AddVertex(nv);
                }
                return nv;
            }

            foreach (Vertex v in Verts.ToArray()) Map(v);

            List<Polygon> polys = new List<Polygon>(Polygons.Count);
            foreach (Polygon p in Polygons)
            {
                Vertex[] src = p.GetVertices();
                Vertex[] dst = new Vertex[src.Length];
                for (int i = 0; i < src.Length; i++) dst[i] = Map(src[i]);
                polys.Add(new Polygon(dst)); // ctor re-links vertex -> polygon associations
            }

            // Copy the cached dual too. Dual polygons own their own Vertex objects (no
            // shared refs with the primal graph), so we build a fresh independent set —
            // consumer mirrors then never see worker-side mutations to vertex positions.
            List<Polygon> dualCopy = null;
            if (Dual != null)
            {
                dualCopy = new List<Polygon>(Dual.Count);
                foreach (Polygon p in Dual)
                {
                    Vertex[] src = p.GetVertices();
                    Vertex[] dst = new Vertex[src.Length];
                    for (int i = 0; i < src.Length; i++)
                        dst[i] = new Vertex(src[i].Position, src[i].IsEdge);
                    dualCopy.Add(new Polygon(dst));
                }
            }

            return new PrimalChunk(Coord, polys, vc)
            {
                Version = Version,
                Dual = dualCopy,
                DualBuiltFromVersion = DualBuiltFromVersion,
            };
        }
    }
}
