using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid
{
    // A regeneration-free snapshot of a chunk's primal grid: vertex positions + edge flags,
    // and polygons as vertex-index tuples. It captures the chunk's *relaxed* state, so a
    // reload restores it exactly rather than regenerating a fresh lattice-aligned chunk
    // (which would no longer match a neighbour that stayed loaded).
    //
    // Stores only primitive arrays — far lighter than the live Vertex/Polygon object graph,
    // and trivially serializable to disk later.
    public class ChunkSnapshot
    {
        public Vector3[] VertexPositions;
        public bool[] VertexIsEdge;
        public int[][] PolygonVertices; // each polygon as indices into the vertex arrays
        public int Version;

        public static ChunkSnapshot Capture(PrimalChunk chunk)
        {
            List<Vertex> indexed = new List<Vertex>();
            Dictionary<Vertex, int> index = new Dictionary<Vertex, int>();

            void Add(Vertex v)
            {
                if (!index.ContainsKey(v)) { index[v] = indexed.Count; indexed.Add(v); }
            }

            // Index every vertex (collection first, then any referenced by polygons).
            foreach (Vertex v in chunk.Verts.ToArray()) Add(v);
            foreach (Polygon p in chunk.Polygons)
                foreach (Vertex v in p.GetVertices()) Add(v);

            ChunkSnapshot snap = new ChunkSnapshot
            {
                VertexPositions = new Vector3[indexed.Count],
                VertexIsEdge = new bool[indexed.Count],
                PolygonVertices = new int[chunk.Polygons.Count][],
                Version = chunk.Version,
            };

            for (int i = 0; i < indexed.Count; i++)
            {
                snap.VertexPositions[i] = indexed[i].Position;
                snap.VertexIsEdge[i] = indexed[i].IsEdge;
            }

            for (int pi = 0; pi < chunk.Polygons.Count; pi++)
            {
                Vertex[] pv = chunk.Polygons[pi].GetVertices();
                int[] ids = new int[pv.Length];
                for (int k = 0; k < pv.Length; k++) ids[k] = index[pv[k]];
                snap.PolygonVertices[pi] = ids;
            }

            return snap;
        }

        public static PrimalChunk Restore(ChunkCoord coord, ChunkSnapshot snap)
        {
            int n = snap.VertexPositions.Length;
            Vertex[] verts = new Vertex[n];
            VertexCollection vc = new VertexCollection();
            for (int i = 0; i < n; i++)
            {
                verts[i] = new Vertex(snap.VertexPositions[i], snap.VertexIsEdge[i]);
                vc.AddVertex(verts[i]);
            }

            List<Polygon> polygons = new List<Polygon>(snap.PolygonVertices.Length);
            foreach (int[] ids in snap.PolygonVertices)
            {
                Vertex[] pv = new Vertex[ids.Length];
                for (int k = 0; k < ids.Length; k++) pv[k] = verts[ids[k]];
                polygons.Add(new Polygon(pv)); // ctor re-links vertex -> polygon associations
            }

            return new PrimalChunk(coord, polygons, vc) { Version = snap.Version };
        }
    }
}
