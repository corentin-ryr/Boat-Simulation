using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    public class PolygonGridGenerator
    {
        // Full pipeline — kept for standalone (non-chunked) use
        public static (HashSet<Polygon>, VertexCollection) GenerateGrid(int gridSize, float hexagonRadius, System.Random random, int nbIterRelaxation, bool normalizedRelaxation)
        {
            VertexCollection primalVerts = new VertexCollection();
            GeneratePrimalPolygons(gridSize, hexagonRadius, random, Vector3.zero, nbIterRelaxation, normalizedRelaxation, primalVerts);
            return GenerateDual(primalVerts);
        }

        // Step 1 — generate this chunk's primal grid (hex → merge → relax), with a world offset.
        // Returns the polygon list and the vertex collection for this chunk only.
        public static (List<Polygon>, VertexCollection) GeneratePrimal(int gridSize, float hexagonRadius, System.Random random, Vector3 offset, int nbIterRelaxation, bool normalizedRelaxation)
        {
            VertexCollection vertices = new VertexCollection();
            List<Polygon> polygons = GeneratePrimalPolygons(gridSize, hexagonRadius, random, offset, nbIterRelaxation, normalizedRelaxation, vertices);
            return (polygons, vertices);
        }

        private static List<Polygon> GeneratePrimalPolygons(int gridSize, float hexagonRadius, System.Random random, Vector3 offset, int nbIterRelaxation, bool normalizedRelaxation, VertexCollection vertices)
        {
            List<Polygon> triangles = GenHexShape(gridSize, hexagonRadius, vertices, offset);
            ComputeNeighbors(triangles);
            triangles = RandomTriangleMerge(random, triangles, vertices);
            ComputeNeighbors(triangles);
            GridRelaxation(nbIterRelaxation, normalizedRelaxation, triangles, vertices);
            return triangles;
        }

        // Generate the dual grid from a chunk's primal vertex collection (standalone / no
        // neighbours — boundary cells will be incomplete; used by GridGenerator).
        public static (HashSet<Polygon>, VertexCollection) GenerateDual(VertexCollection primalVerts)
        {
            (HashSet<Polygon> dual, VertexCollection dualVerts) = GenerateDualGrid(primalVerts);
            ComputeNeighbors(dual.ToList());
            return (dual, dualVerts);
        }

        // Cross-chunk dual: builds this chunk's dual cells.
        //  • Completion — a border (IsEdge) vertex's cell is closed using the neighbour faces
        //    incident to the matching border vertex (looked up by rounded position).
        //  • Ownership — a border vertex shared by several chunks is built by exactly one of
        //    them (lowest ChunkCoord), so seam cells render once with no overlap.
        // `neighborFaces` and `minNeighborCoord` are keyed by LatticeKey(position); the latter
        // holds the lowest coord among the *neighbour* chunks sharing each border position.
        //
        // Allocation-frugal: walks VertexCollection.Values directly, reuses a pair of stack
        // buffers across all vertices (only resized on the rare wide ring), and does the
        // by-angle sort in place via insertion sort. We deliberately do NOT call
        // ComputeNeighbors on the resulting dual polygons — they're consumed by mesh build
        // (no neighbour access) and discarded by DeepCopy (which drops them anyway), so the
        // edge-pairing pass would be pure waste.
        public static (List<Polygon>, VertexCollection) GenerateDual(
            PrimalChunk chunk, float hexRadius,
            Dictionary<(int, int), List<Polygon>> neighborFaces,
            Dictionary<(int, int), ChunkCoord> minNeighborCoord)
        {
            // Dual polygon count is bounded above by primal vertex count.
            List<Polygon> dualPolygons = new List<Polygon>(chunk.Verts.Count);
            VertexCollection dualVertices = new VertexCollection();

            // Reusable buffers. Typical dual cell has 3-7 faces; 16 covers the wide ring
            // around a 3-chunk corner with comfortable headroom.
            Polygon[] facesBuf = new Polygon[16];
            Vertex[]  vertsBuf = new Vertex[16];
            float[]   anglesBuf = new float[16];

            foreach (Vertex vertex in chunk.Verts.Values)
            {
                // Gather incident faces into the local buffer (no per-vertex List allocation).
                int faceCount = 0;
                foreach (Polygon p in vertex.Polygons)
                {
                    if (faceCount >= facesBuf.Length) Array.Resize(ref facesBuf, facesBuf.Length * 2);
                    facesBuf[faceCount++] = p;
                }

                if (vertex.IsEdge)
                {
                    (int, int) key = LatticeKey(vertex.Position, hexRadius);

                    // Ownership: skip if a sharing neighbour has a lower coord — it builds it.
                    if (minNeighborCoord != null
                        && minNeighborCoord.TryGetValue(key, out ChunkCoord minN)
                        && minN.CompareTo(chunk.Coord) < 0)
                        continue;

                    // Completion: add the neighbour faces so the ring closes across the seam.
                    if (neighborFaces != null && neighborFaces.TryGetValue(key, out List<Polygon> extra))
                    {
                        for (int i = 0, n = extra.Count; i < n; i++)
                        {
                            if (faceCount >= facesBuf.Length) Array.Resize(ref facesBuf, facesBuf.Length * 2);
                            facesBuf[faceCount++] = extra[i];
                        }
                    }
                }

                if (faceCount < 3) continue; // exposed/degenerate — drop instead of spiking

                if (faceCount > vertsBuf.Length)
                {
                    Array.Resize(ref vertsBuf, faceCount);
                    Array.Resize(ref anglesBuf, faceCount);
                }

                // Materialize dual vertices (face centers) and their sort key (angle from
                // this primal vertex around the up axis) in one pass.
                Vector3 vp = vertex.Position;
                for (int i = 0; i < faceCount; i++)
                {
                    Vertex dv = dualVertices.AddOrCreate(facesBuf[i].GetCenter());
                    vertsBuf[i] = dv;
                    anglesBuf[i] = Vector3.SignedAngle(dv.Position - vp, Vector3.forward, Vector3.up);
                }

                // Insertion sort by angle — faceCount is small (3-8 normally, up to ~12
                // around a 3-chunk corner) so the simpler algorithm beats a comparator-based
                // Array.Sort that would allocate a delegate per call.
                for (int i = 1; i < faceCount; i++)
                {
                    float a = anglesBuf[i];
                    Vertex v = vertsBuf[i];
                    int j = i - 1;
                    while (j >= 0 && anglesBuf[j] > a)
                    {
                        anglesBuf[j + 1] = anglesBuf[j];
                        vertsBuf[j + 1] = vertsBuf[j];
                        j--;
                    }
                    anglesBuf[j + 1] = a;
                    vertsBuf[j + 1] = v;
                }

                // Polygon takes ownership of the vertex array, so we must allocate a
                // right-sized copy rather than handing it the shared buffer.
                Vertex[] polyVerts = new Vertex[faceCount];
                Array.Copy(vertsBuf, polyVerts, faceCount);
                dualPolygons.Add(new Polygon(polyVerts));
            }

            return (dualPolygons, dualVertices);
        }

        #region Border relaxation

        // Joint relaxation of shared border (IsEdge) vertices across loaded chunks, plus
        // `interiorDepth` rings of interior vertices just inside each active seam.
        //
        // Seam vertices that exist in two or more chunks are moved using the combined quad
        // neighbourhood from every chunk that shares them, and the same delta is applied to
        // each copy so the seam stays welded (2-chunk edges and 3-chunk corners are handled
        // uniformly — a corner just has three copies). Relaxing only the seam leaves a kink
        // at the first interior ring (it was relaxed while the border was pinned), so the
        // next `interiorDepth` rings of non-shared interior vertices are also relaxed, each
        // using its own chunk's local neighbourhood. Vertices deeper than that, and exposed
        // boundaries with no loaded neighbour, stay put.
        // Returns the indices (into `chunks`) of every chunk whose vertices actually moved
        // beyond a small epsilon, so the caller can rebuild only those chunks' meshes.
        public static HashSet<int> RelaxBorders(IReadOnlyList<VertexCollection> chunks, float hexRadius, int iterations, int interiorDepth = 1)
        {
            HashSet<int> movedChunks = new HashSet<int>();
            const float epsSq = 1e-6f; // (1e-3)^2 — sub-millimetre at unit scale

            // Group every border vertex copy by its deterministic lattice index, and remember
            // which chunk each copy belongs to (owner) for version attribution. Keying on the
            // exact integer index (not a rounded position) guarantees coincident copies in
            // adjacent chunks land in the same group regardless of float error.
            Dictionary<Vertex, int> owner = new Dictionary<Vertex, int>();
            Dictionary<(int, int), List<Vertex>> groups = new Dictionary<(int, int), List<Vertex>>();
            for (int ci = 0; ci < chunks.Count; ci++)
                foreach (Vertex v in chunks[ci].ToArray())
                {
                    if (!v.IsEdge) continue;
                    owner[v] = ci;
                    (int, int) key = LatticeKey(v.Position, hexRadius);
                    if (!groups.TryGetValue(key, out List<Vertex> list))
                        groups[key] = list = new List<Vertex>();
                    list.Add(v);
                }

            // Seam groups (moved jointly). Only positions shared by 2+ loaded chunks.
            List<List<Vertex>> seams = groups.Values.Where(g => g.Count >= 2).ToList();
            if (seams.Count == 0) return movedChunks;

            // BFS outward from the seam to collect interiorDepth rings of interior vertices
            // (moved individually, never shared). Each interior vertex inherits the owner of
            // the seam/interior vertex it was reached from (same chunk).
            HashSet<Vertex> visited = new HashSet<Vertex>();
            foreach (List<Vertex> g in seams)
                foreach (Vertex v in g) visited.Add(v);

            List<Vertex> interior = new List<Vertex>();
            List<Vertex> frontier = visited.ToList();
            for (int d = 0; d < interiorDepth; d++)
            {
                List<Vertex> next = new List<Vertex>();
                foreach (Vertex v in frontier)
                {
                    int ownerCi = owner.TryGetValue(v, out int o) ? o : -1;
                    foreach (Polygon p in v.Polygons)
                        foreach (Vertex u in p.GetVertices())
                            if (!u.IsEdge && visited.Add(u))
                            {
                                if (ownerCi >= 0) owner[u] = ownerCi;
                                interior.Add(u);
                                next.Add(u);
                            }
                }
                frontier = next;
            }

            for (int it = 0; it < iterations; it++)
            {
                // Compute every delta before applying so each iteration is simultaneous.
                Vector3[] seamDeltas = new Vector3[seams.Count];
                for (int i = 0; i < seams.Count; i++)
                {
                    Vector3 sum = Vector3.zero;
                    foreach (Vertex v in seams[i]) sum += VertexGoalMovement(v);
                    seamDeltas[i] = sum * 0.1f;
                }

                Vector3[] interiorDeltas = new Vector3[interior.Count];
                for (int i = 0; i < interior.Count; i++)
                    interiorDeltas[i] = VertexGoalMovement(interior[i]) * 0.1f;

                // Apply a delta only if it clears the epsilon, so translation and the
                // version bump are atomic: a seam either moves >eps (and every copy moves by
                // the identical delta AND every owning chunk is marked, so all copies are
                // re-published together) or it doesn't move at all. Translating sub-eps would
                // let two copies of a seam be snapshotted at different states across passes
                // (only the chunk whose version happened to bump gets re-captured), drifting
                // them apart until the cross-chunk dual lookup misses — the grey seam ribbons.
                for (int i = 0; i < seams.Count; i++)
                {
                    if (seamDeltas[i].sqrMagnitude <= epsSq) continue;
                    foreach (Vertex v in seams[i])
                    {
                        v.Translate(seamDeltas[i]);
                        if (owner.TryGetValue(v, out int ci)) movedChunks.Add(ci);
                    }
                }
                for (int i = 0; i < interior.Count; i++)
                {
                    if (interiorDeltas[i].sqrMagnitude <= epsSq) continue;
                    interior[i].Translate(interiorDeltas[i]);
                    if (owner.TryGetValue(interior[i], out int ci)) movedChunks.Add(ci);
                }
            }

            return movedChunks;
        }

        // Rotational (square-fitting) relaxation goal contributed to vertex v by every quad
        // incident to it — the same formula as GridRelaxation, summed for a single vertex.
        //
        // Returns XZ-only: relaxation reshapes the planar layout of the grid, while Y carries
        // elevation sampled from the noise field. The rotations are around Vector3.up so they
        // preserve Y on each rotated component, but the SUM still inherits Y from the input
        // vertex positions; zeroing it here keeps relaxation conceptually 2D and elevation
        // untouched by seam welding.
        private static Vector3 VertexGoalMovement(Vertex v)
        {
            Vector3 sum = Vector3.zero;
            foreach (Polygon p in v.Polygons)
            {
                Vertex[] verts = p.GetVertices();
                if (verts.Length != 4) continue;

                Vector3 c = p.GetCenter();
                Vector3 goal = (verts[0].Position - c
                              + Quaternion.AngleAxis(-270, Vector3.up) * (verts[1].Position - c)
                              + Quaternion.AngleAxis(-180, Vector3.up) * (verts[2].Position - c)
                              + Quaternion.AngleAxis(-90, Vector3.up) * (verts[3].Position - c)) / 4f;

                int i = System.Array.IndexOf(verts, v);
                Vector3 target = c + Quaternion.AngleAxis(-90f * i, Vector3.up) * goal;
                sum += target - v.Position;
            }
            sum.y = 0f;
            return sum;
        }

        // Deterministic lattice identity for a border vertex, used to match the same seam point
        // across adjacent chunks. Border vertices lie on the globally-shared triangle half-lattice
        // (lattice points + ortho-subdivision edge midpoints): x = (I + J/2)·s, z = J·(√3/2)·s
        // with I,J on the half-integer grid. Returning the doubled indices (2I, 2J) as integers
        // makes coincident copies in two chunks key *identically* — exact, with none of a rounded
        // position's bucket-boundary fragility: nodes are 0.5·s apart and a vertex sits within a
        // few ULPs of its node, so the rounding is never close to a decision boundary.
        public static (int, int) LatticeKey(Vector3 p, float hexRadius)
        {
            float h = Mathf.Sqrt(3f) / 2f;
            int j2 = Mathf.RoundToInt(2f * p.z / (h * hexRadius));
            int i2 = Mathf.RoundToInt(2f * p.x / hexRadius - j2 * 0.5f);
            return (i2, j2);
        }

        #endregion

        #region Hexagon grid creation

        // Generates a hexagonal region of a single uniform equilateral-triangle lattice
        // (Townscaper-style base grid). Triangle side length = hexagonRadius. All vertices
        // sit on one triangular lattice — there are no per-cell center vertices. The outer
        // boundary is flagged later by ComputeNeighbors.
        private static List<Polygon> GenHexShape(int gridSize, float hexagonRadius, VertexCollection vertices, Vector3 offset = default)
        {
            List<Polygon> triangles = new List<Polygon>();
            float s = hexagonRadius;
            float h = Mathf.Sqrt(3f) / 2f;

            bool InRegion(int i, int j) =>
                Mathf.Abs(i) <= gridSize && Mathf.Abs(j) <= gridSize && Mathf.Abs(i + j) <= gridSize;

            Vertex Lattice(int i, int j)
            {
                float x = (i + j * 0.5f) * s + offset.x;
                float z = j * h * s + offset.z;
                return vertices.AddOrCreate(new Vector3(x, 0f, z));
            }

            for (int i = -gridSize; i <= gridSize; i++)
            {
                for (int j = -gridSize; j <= gridSize; j++)
                {
                    // Up triangle: (i,j), (i+1,j), (i,j+1)
                    if (InRegion(i, j) && InRegion(i + 1, j) && InRegion(i, j + 1))
                        triangles.Add(new Polygon(new Vertex[]
                        {
                            Lattice(i, j), Lattice(i + 1, j), Lattice(i, j + 1)
                        }));

                    // Down triangle: (i+1,j), (i+1,j+1), (i,j+1)
                    if (InRegion(i + 1, j) && InRegion(i + 1, j + 1) && InRegion(i, j + 1))
                        triangles.Add(new Polygon(new Vertex[]
                        {
                            Lattice(i + 1, j), Lattice(i + 1, j + 1), Lattice(i, j + 1)
                        }));
                }
            }

            return triangles;
        }

        // Build polygon adjacency by indexing each undirected edge once: every edge is shared by
        // at most two polygons, so a single pass over all edges pairs them up. O(total edges)
        // with no per-vertex/per-polygon temporary allocations (the old version was a quadratic
        // candidate scan that re-allocated arrays on every vertex access). Edges left unpaired
        // are the boundary, and their vertices are flagged IsEdge.
        private static void ComputeNeighbors(List<Polygon> polygons)
        {
            foreach (Polygon polygon in polygons) polygon.InitNeighbors();

            // edge -> the first polygon seen on it (and which edge index). Matched edges are
            // removed, so whatever remains at the end is boundary.
            var pending = new Dictionary<(Vertex, Vertex), (Polygon poly, int edge)>();

            foreach (Polygon polygon in polygons)
            {
                Vertex[] verts = polygon.GetVertices();
                for (int i = 0; i < verts.Length; i++)
                {
                    (Vertex, Vertex) key = EdgeKey(verts[i], verts[(i + 1) % verts.Length]);
                    if (pending.TryGetValue(key, out (Polygon poly, int edge) other))
                    {
                        polygon.SetNeighborAt(i, other.poly);
                        other.poly.SetNeighborAt(other.edge, polygon);
                        pending.Remove(key);
                    }
                    else pending[key] = (polygon, i);
                }
            }

            foreach (KeyValuePair<(Vertex, Vertex), (Polygon, int)> kv in pending)
            {
                kv.Key.Item1.IsEdge = true;
                kv.Key.Item2.IsEdge = true;
            }
        }

        // Order-independent key for an undirected edge: order the two endpoints by position so
        // both incident polygons (whatever their winding) produce the same key. Distinct welded
        // vertices have distinct positions, so this is unambiguous.
        private static (Vertex, Vertex) EdgeKey(Vertex a, Vertex b)
        {
            Vector3 pa = a.Position, pb = b.Position;
            bool aFirst = pa.x < pb.x || (pa.x == pb.x && pa.z < pb.z);
            return aFirst ? (a, b) : (b, a);
        }

        #endregion

        #region Convert hexagons to quads

        private static List<Polygon> RandomTriangleMerge(System.Random random, List<Polygon> triangles, VertexCollection vertices)
        {
            List<Polygon> polygonsAfterMerge = new List<Polygon>();

            while (triangles.Count > 0)
            {
                int index = random.Next(triangles.Count);
                Polygon current = triangles[index];

                // Border triangles ARE allowed to merge with interior neighbours, so the
                // randomness reaches the seams instead of leaving a regular band. Boundary
                // edges have no neighbour and never merge, so the seam vertices (lattice
                // corners + their edge midpoints) stay deterministic and coincide across
                // chunks regardless of how interior merging resolves on each side.
                List<Polygon> triangleNeighbors = current.GetNeighbors()
                    .Where(n => n != null && n.GetVertices().Length == 3)
                    .ToList();

                if (triangleNeighbors.Count > 0)
                {
                    Polygon neighbor = triangleNeighbors[random.Next(triangleNeighbors.Count)];
                    Vertex[] tempVertices = current.GetVertices();

                    HashSet<Vertex> union = tempVertices.ToHashSet();
                    union.UnionWith(neighbor.GetVertices());

                    HashSet<Vertex> intersect = tempVertices.ToHashSet();
                    intersect.IntersectWith(neighbor.GetVertices());
                    union.ExceptWith(intersect);

                    Vertex[] outer = union.ToArray();

                    Vertex first, second;
                    if (outer[0] == tempVertices[0]) { first = tempVertices[1]; second = tempVertices[2]; }
                    else if (outer[0] == tempVertices[1]) { first = tempVertices[2]; second = tempVertices[0]; }
                    else { first = tempVertices[0]; second = tempVertices[1]; }

                    Polygon square = new Polygon(new Vertex[] { outer[0], first, outer[1], second });

                    foreach (Vertex v in union) v.AddPolygon(square);
                    foreach (Vertex v in intersect) v.AddPolygon(square);

                    foreach (Polygon n in current.GetNeighbors()) n?.UpdateNeigbhor(current, null);
                    foreach (Polygon n in neighbor.GetNeighbors()) n?.UpdateNeigbhor(neighbor, null);

                    polygonsAfterMerge.AddRange(SubdividePolygon(square, vertices));

                    foreach (Vertex v in tempVertices) v.RemovePolygon(current);
                    foreach (Vertex v in neighbor.GetVertices()) v.RemovePolygon(neighbor);

                    triangles.Remove(neighbor);
                }
                else
                {
                    polygonsAfterMerge.AddRange(SubdividePolygon(current, vertices));
                }

                triangles.Remove(current);
            }

            return polygonsAfterMerge;
        }

        private static List<Polygon> SubdividePolygon(Polygon polygon, VertexCollection vertices)
        {
            Vertex center = vertices.AddOrCreate(polygon.GetCenter());
            List<Polygon> subdivided = new List<Polygon>();

            Vertex[] verts = polygon.GetVertices();
            List<Vertex> edgeMidpoints = new List<Vertex>();
            for (int i = 0; i < verts.Length; i++)
                edgeMidpoints.Add(vertices.AddOrCreate((verts[i].Position + verts[(i + 1) % verts.Length].Position) / 2));

            for (int i = 0; i < verts.Length; i++)
            {
                int prev = (i - 1 + verts.Length) % verts.Length;
                subdivided.Add(new Polygon(new Vertex[] { verts[i], edgeMidpoints[i], center, edgeMidpoints[prev] }));
                verts[i].RemovePolygon(polygon);
            }

            return subdivided;
        }

        private static void GridRelaxation(int nbIterRelaxation, bool normalizedRelaxation, List<Polygon> polygons, VertexCollection vertices)
        {
            for (int i = 0; i < nbIterRelaxation; i++)
            {
                float avgDist = 0f;
                foreach (Polygon p in polygons)
                {
                    Vector3 c = p.GetCenter();
                    foreach (Vertex v in p.GetVertices()) avgDist += (c - v.Position).magnitude;
                }
                avgDist /= 4 * polygons.Count;

                foreach (Polygon p in polygons)
                {
                    Vertex[] verts = p.GetVertices();
                    Vector3 c = p.GetCenter();

                    Vector3 goal = (verts[0].Position - c +
                                   Quaternion.AngleAxis(-270, Vector3.up) * (verts[1].Position - c) +
                                   Quaternion.AngleAxis(-180, Vector3.up) * (verts[2].Position - c) +
                                   Quaternion.AngleAxis(-90, Vector3.up) * (verts[3].Position - c)) / 4;

                    if (normalizedRelaxation) goal /= goal.magnitude * avgDist;

                    verts[0].AccumulateMovement(c + goal - verts[0].Position);
                    verts[1].AccumulateMovement(c + Quaternion.AngleAxis(270, Vector3.up) * goal - verts[1].Position);
                    verts[2].AccumulateMovement(c + Quaternion.AngleAxis(180, Vector3.up) * goal - verts[2].Position);
                    verts[3].AccumulateMovement(c + Quaternion.AngleAxis(90, Vector3.up) * goal - verts[3].Position);
                }

                foreach (Vertex v in vertices.ToArray()) v.UpdateVertexPosition();
            }
        }

        private static (HashSet<Polygon>, VertexCollection) GenerateDualGrid(VertexCollection vertices)
        {
            HashSet<Polygon> dualPolygons = new HashSet<Polygon>();
            VertexCollection dualVertices = new VertexCollection();

            foreach (Vertex vertex in vertices.ToArray())
            {
                List<Vertex> newVerts = vertex.Polygons
                    .Select(p => dualVertices.AddOrCreate(p.GetCenter()))
                    .ToList();

                if (newVerts.Count == 2) newVerts.Add(dualVertices.AddOrCreate(vertex.Position));

                newVerts = newVerts
                    .OrderBy(x => Vector3.SignedAngle(x.Position - vertex.Position, Vector3.forward, Vector3.up))
                    .ToList();

                dualPolygons.Add(new Polygon(newVerts.ToArray()));
            }

            return (dualPolygons, dualVertices);
        }

        #endregion
    }
}
