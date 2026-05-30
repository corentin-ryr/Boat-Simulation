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
        public Vector3[] VertexPositions;       // includes Y from the elevation field
        public bool[] VertexIsEdge;
        public int[][] PolygonVertices;         // each polygon as indices into the vertex arrays
        public CellTerrain[] PolygonTerrain;    // per-polygon classification (Ocean/Coastal/Land)
        public bool IsFlat;                     // chunk-level flat-ocean fast-path flag
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
                PolygonTerrain = new CellTerrain[chunk.Polygons.Count],
                IsFlat = chunk.IsFlat,
                Version = chunk.Version,
            };

            for (int i = 0; i < indexed.Count; i++)
            {
                snap.VertexPositions[i] = indexed[i].Position;
                snap.VertexIsEdge[i] = indexed[i].IsEdge;
            }

            for (int pi = 0; pi < chunk.Polygons.Count; pi++)
            {
                Polygon poly = chunk.Polygons[pi];
                Vertex[] pv = poly.GetVertices();
                int[] ids = new int[pv.Length];
                for (int k = 0; k < pv.Length; k++) ids[k] = index[pv[k]];
                snap.PolygonVertices[pi] = ids;
                snap.PolygonTerrain[pi] = poly.Terrain;
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
            for (int pi = 0; pi < snap.PolygonVertices.Length; pi++)
            {
                int[] ids = snap.PolygonVertices[pi];
                Vertex[] pv = new Vertex[ids.Length];
                for (int k = 0; k < ids.Length; k++) pv[k] = verts[ids[k]];
                Polygon poly = new Polygon(pv); // ctor re-links vertex -> polygon associations
                // Older snapshots predating the elevation system have no Terrain array; default
                // to Ocean (consistent with the pre-elevation flat-Y=0 state).
                poly.Terrain = snap.PolygonTerrain != null && pi < snap.PolygonTerrain.Length
                    ? snap.PolygonTerrain[pi]
                    : CellTerrain.Ocean;
                polygons.Add(poly);
            }

            return new PrimalChunk(coord, polygons, vc)
            {
                Version = snap.Version,
                IsFlat = snap.IsFlat,
            };
        }
    }
}
