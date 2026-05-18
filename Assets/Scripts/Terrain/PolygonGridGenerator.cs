using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    public class PolygonGridGenerator
    {
        public static (HashSet<Polygon>, VertexCollection) GenerateGrid(int gridSize, float hexagonRadius, System.Random random, int nbIterRelaxation, bool normalizedRelaxation)
        {
            VertexCollection vertices = new VertexCollection();

            List<Polygon> triangles = GenHexShape(gridSize, hexagonRadius, vertices);
            ComputeNeighbors(triangles);
            triangles = RandomTriangleMerge(random, triangles, vertices);
            ComputeNeighbors(triangles);
            GridRelaxation(nbIterRelaxation, normalizedRelaxation, triangles, vertices);

            (HashSet<Polygon> dualPolygons, VertexCollection dualVertices) = GenerateDualGrid(vertices);
            ComputeNeighbors(dualPolygons.ToList());

            return (dualPolygons, dualVertices);
        }

        #region Hexagon grid creation

        private static List<Polygon> GenHexShape(int gridSize, float hexagonRadius, VertexCollection vertices)
        {
            List<Polygon> triangles = new List<Polygon>();

            for (int q = -gridSize; q <= gridSize; q++)
            {
                int r1 = Mathf.Max(-gridSize, -q - gridSize);
                int r2 = Mathf.Min(gridSize, -q + gridSize);
                for (int r = r1; r <= r2; r++)
                {
                    double height = 2 * hexagonRadius;
                    double width = Math.Sqrt(3) * hexagonRadius;

                    double posx = hexagonRadius * Math.Sqrt(3.0) * (q + r / 2.0);
                    double posz = hexagonRadius * 3.0 / 2.0 * r;

                    Vertex centerVertex = vertices.AddOrCreate(new Vector3((float)posx, 0, (float)posz));
                    Vertex vertex1 = vertices.AddOrCreate(new Vector3((float)posx, 0, (float)(posz + height / 2.0)));
                    Vertex vertex2 = vertices.AddOrCreate(new Vector3((float)(posx + width / 2.0), 0, (float)(posz + height / 4.0)));
                    Vertex vertex3 = vertices.AddOrCreate(new Vector3((float)(posx + width / 2.0), 0, (float)(posz - height / 4.0)));
                    Vertex vertex4 = vertices.AddOrCreate(new Vector3((float)posx, 0, (float)(posz - height / 2.0)));
                    Vertex vertex5 = vertices.AddOrCreate(new Vector3((float)(posx - width / 2.0), 0, (float)(posz - height / 4.0)));
                    Vertex vertex6 = vertices.AddOrCreate(new Vector3((float)(posx - width / 2.0), 0, (float)(posz + height / 4.0)));

                    triangles.AddRange(new Polygon[] {
                        new Polygon(new Vertex[] { vertex1, centerVertex, vertex2 }),
                        new Polygon(new Vertex[] { vertex2, centerVertex, vertex3 }),
                        new Polygon(new Vertex[] { vertex3, centerVertex, vertex4 }),
                        new Polygon(new Vertex[] { vertex4, centerVertex, vertex5 }),
                        new Polygon(new Vertex[] { vertex5, centerVertex, vertex6 }),
                        new Polygon(new Vertex[] { vertex6, centerVertex, vertex1 }),
                    });
                }
            }

            return triangles;
        }

        private static void ComputeNeighbors(List<Polygon> polygons)
        {
            foreach (Polygon polygon in polygons)
            {
                HashSet<Polygon> neighborCandidates = new HashSet<Polygon>();
                foreach (Vertex vertex in polygon.GetVertices()) neighborCandidates.UnionWith(vertex.Polygons);
                neighborCandidates.Remove(polygon);

                List<Polygon> neighbors = new List<Polygon>();
                Vertex[] verts = polygon.GetVertices();
                for (int i = 0; i < verts.Length; i++)
                {
                    Vertex v1 = verts[i];
                    Vertex v2 = verts[(i + 1) % verts.Length];

                    bool foundNeighbor = false;
                    foreach (Polygon candidate in neighborCandidates)
                    {
                        if (candidate.GetVertices().Contains(v1) && candidate.GetVertices().Contains(v2))
                        {
                            neighbors.Add(candidate);
                            foundNeighbor = true;
                            break;
                        }
                    }
                    if (!foundNeighbor)
                    {
                        neighbors.Add(null);
                        v1.IsEdge = true;
                        v2.IsEdge = true;
                    }
                }
                polygon.SetNeighbors(neighbors.ToArray());
            }
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
