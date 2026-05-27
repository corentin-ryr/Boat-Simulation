using System.Collections.Generic;

namespace TerrainGrid
{
    // One chunk's primal (data-layer) grid. Pure data — no rendering, no Unity objects.
    public class PrimalChunk
    {
        public readonly ChunkCoord Coord;
        public List<Polygon> Polygons;
        public VertexCollection Verts;

        // Incremented by the border-relaxation pass whenever this chunk's vertices move,
        // so the render layer can detect a stale mesh and rebuild only when needed.
        public int Version;

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

            return new PrimalChunk(Coord, polys, vc) { Version = Version };
        }
    }
}
