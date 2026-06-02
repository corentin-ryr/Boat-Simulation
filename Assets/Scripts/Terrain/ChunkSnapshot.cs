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
        public TileProperty[] PolygonProperties; // per-polygon gameplay slot; null = all Default
        public bool IsFlat;                     // chunk-level flat-ocean fast-path flag
        public int Version;

        // Cached dual, persisted so a restore behaves like a cache revive — without these
        // fields, ChunkSnapshot would drop a chunk's already-built Dual on the floor and the
        // next pass would re-run BuildDual to produce the same result. Dual polygons own their
        // own Vertex objects (no sharing across polygons — see PrimalChunk.DeepCopy), so we
        // store positions and edge flags inline per polygon instead of indexing into a flat
        // vertex array. Null when the chunk hadn't built its dual at capture time (e.g. fresh
        // gen evicted before first pass) OR for flat chunks (their dual is the empty sentinel,
        // reconstructable from IsFlat alone).
        public Vector3[][] DualPolygonPositions;
        public bool[][] DualPolygonIsEdge;
        public bool DualComplete;
        public int DualBuiltFromVersion = -1;

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

            // Persist gameplay tile slots only when actually populated; null in the snap
            // means "every cell defaulted" and restores to a null array on the chunk too.
            if (chunk.TileProperties != null)
            {
                int n = chunk.Polygons.Count;
                snap.PolygonProperties = new TileProperty[n];
                System.Array.Copy(chunk.TileProperties, snap.PolygonProperties, n);
            }

            // Persist the cached dual for non-flat chunks. Flat chunks' dual is the shared
            // empty sentinel — reconstructable from IsFlat alone, no need to write it out.
            if (chunk.Dual != null && !chunk.IsFlat)
            {
                int dCount = chunk.Dual.Count;
                snap.DualPolygonPositions = new Vector3[dCount][];
                snap.DualPolygonIsEdge = new bool[dCount][];
                for (int pi = 0; pi < dCount; pi++)
                {
                    Vertex[] dv = chunk.Dual[pi].GetVertices();
                    Vector3[] pos = new Vector3[dv.Length];
                    bool[] edge = new bool[dv.Length];
                    for (int k = 0; k < dv.Length; k++)
                    {
                        pos[k] = dv[k].Position;
                        edge[k] = dv[k].IsEdge;
                    }
                    snap.DualPolygonPositions[pi] = pos;
                    snap.DualPolygonIsEdge[pi] = edge;
                }
                snap.DualComplete = chunk.DualComplete;
                snap.DualBuiltFromVersion = chunk.DualBuiltFromVersion;
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

            PrimalChunk pc = new PrimalChunk(coord, polygons, vc)
            {
                Version = snap.Version,
                IsFlat = snap.IsFlat,
            };

            // Gameplay tile slots: restore only when present and length-matched, matching the
            // null-fallback convention used for older snapshots' missing PolygonTerrain above.
            if (snap.PolygonProperties != null && snap.PolygonProperties.Length == polygons.Count)
            {
                pc.TileProperties = new TileProperty[polygons.Count];
                System.Array.Copy(snap.PolygonProperties, pc.TileProperties, polygons.Count);
            }

            if (snap.IsFlat)
            {
                // Flat chunks need no neighbour info — their "dual" is the shared empty
                // sentinel that BuildDual would assign. NeedsDual already short-circuits
                // on IsFlat so we don't bother setting Dual here; but we must mark
                // DualComplete=true so the add-cascade skips this chunk.
                pc.DualComplete = true;
            }
            else if (snap.DualPolygonPositions != null)
            {
                // Round-trip the cached dual: rebuild independent Vertex objects per
                // polygon (matching PrimalChunk.DeepCopy's no-sharing convention). After
                // this, NeedsDual returns false and the chunk skips BuildDual entirely —
                // which is the whole point of persisting the dual.
                int dCount = snap.DualPolygonPositions.Length;
                List<Polygon> dual = new List<Polygon>(dCount);
                bool[][] edges = snap.DualPolygonIsEdge;
                for (int pi = 0; pi < dCount; pi++)
                {
                    Vector3[] pos = snap.DualPolygonPositions[pi];
                    bool[] edge = edges != null && pi < edges.Length ? edges[pi] : null;
                    Vertex[] dv = new Vertex[pos.Length];
                    for (int k = 0; k < pos.Length; k++)
                        dv[k] = new Vertex(pos[k], edge != null && k < edge.Length && edge[k]);
                    dual.Add(new Polygon(dv));
                }
                pc.Dual = dual;
                pc.DualBuiltFromVersion = snap.DualBuiltFromVersion;
                pc.DualComplete = snap.DualComplete;
            }
            // else: older snapshot or chunk that hadn't built its dual yet at capture time.
            // Falls back to the original behaviour (Dual=null → BuildDual will run next pass).

            return pc;
        }
    }
}
