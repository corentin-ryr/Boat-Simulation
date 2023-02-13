using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Habrador_Computational_Geometry;
using UnityEngine;

public class Triangle
{
    public Triangle(int triangleIndex, Vector3 v1, Vector3 v2, Vector3 v3) //Vector3 vertex1, Vector3 vertex2, Vector3 vertex3, 
    {
        this.TriangleIndex = triangleIndex;
        this.Vertex1 = v1;
        this.Vertex2 = v2;
        this.Vertex3 = v3;
    }
    public Triangle(int triangleIndex, Vector3 v1, Vector3 v2, Vector3 v3, int n1, int n2, int n3) //Vector3 vertex1, Vector3 vertex2, Vector3 vertex3, 
    {
        this.TriangleIndex = triangleIndex;
        this.Vertex1 = v1;
        this.Vertex2 = v2;
        this.Vertex3 = v3;

        this.N1 = n1;
        this.N2 = n2;
        this.N3 = n3;
    }

    public int TriangleIndex { get; }

    private int n1;
    private int n2;
    private int n3;
    public int N1 { get => n1; set => n1 = value; }
    public int N2 { get => n2; set => n2 = value; }
    public int N3 { get => n3; set => n3 = value; }

    //The 3 vertices of the triangle
    private Vector3 vertex1;
    private Vector3 vertex2;
    private Vector3 vertex3;

    public Vector3 Vertex1 { get => vertex1; set => vertex1 = value; }
    public Vector3 Vertex2 { get => vertex2; set => vertex2 = value; }
    public Vector3 Vertex3 { get => vertex3; set => vertex3 = value; }
    public void SetVertex1Height(float height) { this.vertex1.y = height; }
    public void SetVertex2Height(float height) { this.vertex2.y = height; }
    public void SetVertex3Height(float height) { this.vertex3.y = height; }
    public Vector3[] GetVertices()
    {
        return new Vector3[] { Vertex1, Vertex2, Vertex3 };
    }


    //The 3 neighbors
    public Triangle T1 { get; set; }//Neighbor between triangle[TriangleIndex * 3] and triangle[TriangleIndex * 3+1]
    public Triangle T2 { get; set; }//Neighbor between triangle[TriangleIndex * 3+1] and triangle[TriangleIndex * 3+2]
    public Triangle T3 { get; set; }//Neighbor between triangle[TriangleIndex * 3+2] and triangle[TriangleIndex * 3]

    public Triangle[] GetNeighbors()
    {
        return new Triangle[] { T1, T2, T3 };
    }

    public override string ToString()
    {
        return "Triangle: " + this.N1 + ", "  + this.N2 + ", "  + this.N3;
    }

}

public struct Cell
{
    public Cell(Bounds bounds, int setCapacity)
    {
        this.Bounds = bounds;
        triangleSet1 = new List<Triangle>(setCapacity);
        triangleSet2 = new List<Triangle>(setCapacity);
        hasCandidatePotential = false;
    }

    public bool Intersects(Bounds bounds)
    {
        return this.Bounds.Intersects(bounds);
    }

    public Bounds Bounds { get; }
    private List<Triangle> triangleSet1;
    public Triangle[] TriangleSet1 { get => triangleSet1.ToArray(); }
    private bool hasCandidatePotential;
    public bool HasCandidatePotential { get => hasCandidatePotential; }

    public void AddSet1(Triangle triangle)
    {
        triangleSet1.Add(triangle);
        hasCandidatePotential = true;
    }
    private List<Triangle> triangleSet2;
    public Triangle[] TriangleSet2 { get => triangleSet2.ToArray(); }

    public void AddSet2(Triangle triangle)
    {
        triangleSet2.Add(triangle);
    }

    public bool HasCandidates()
    {
        if (hasCandidatePotential && triangleSet2.Count > 0)
        {
            return true;
        }
        return false;
    }

    public void resetSet2()
    {
        triangleSet2.Clear();
    }
}

public struct MovingAverage
{
    private List<Vector3> lastFewPoints;
    private int MANumber;
    public MovingAverage(int MANumber)
    {
        lastFewPoints = new List<Vector3>(MANumber);
        this.MANumber = MANumber;
    }

    public Vector3 GetSmoothedValue(Vector3 newValue)
    {
        if (lastFewPoints.Count < MANumber)
        {
            lastFewPoints.Add(newValue);
        }
        else
        {
            lastFewPoints.RemoveAt(0);
            lastFewPoints.Add(newValue);
        }

        return new Vector3(
            lastFewPoints.Average(x => x.x),
            lastFewPoints.Average(x => x.y),
            lastFewPoints.Average(x => x.z)
        );
    }
}

public static class MeshHelper
{

    public static Triangle[] FindTriangleNeighbors(Mesh mesh)
    {
        int Nt = mesh.triangles.Length / 3;
        int Nv = mesh.vertexCount;

        List<Triangle>[] S = new List<Triangle>[Nv];
        Triangle[] triangleNeighbors = new Triangle[Nt];
        for (int i = 0; i < Nt; i++)
        {
            Triangle triangle = new Triangle(i, mesh.vertices[mesh.triangles[i * 3]],
                                                mesh.vertices[mesh.triangles[i * 3 + 1]],
                                                mesh.vertices[mesh.triangles[i * 3 + 2]],
                                                i * 3, i * 3 + 1, i * 3 + 2);
            triangleNeighbors[i] = triangle;
            (S[mesh.triangles[i * 3]] ??= new List<Triangle>()).Add(triangle);
            (S[mesh.triangles[i * 3 + 1]] ??= new List<Triangle>()).Add(triangle);
            (S[mesh.triangles[i * 3 + 2]] ??= new List<Triangle>()).Add(triangle);
        }

        for (int i = 0; i < Nt; i++)
        {
            IEnumerable<Triangle> t1 = S[mesh.triangles[i * 3]].Intersect(S[mesh.triangles[i * 3 + 1]]);
            IEnumerable<Triangle> t2 = S[mesh.triangles[i * 3 + 1]].Intersect(S[mesh.triangles[i * 3 + 2]]);
            IEnumerable<Triangle> t3 = S[mesh.triangles[i * 3 + 2]].Intersect(S[mesh.triangles[i * 3]]);

            triangleNeighbors[i].T1 = t1.Any(p => p.TriangleIndex != i) ? t1.First<Triangle>(p => p.TriangleIndex != i) : null;
            triangleNeighbors[i].T2 = t2.Any(p => p.TriangleIndex != i) ? t2.First<Triangle>(p => p.TriangleIndex != i) : null;
            triangleNeighbors[i].T3 = t3.Any(p => p.TriangleIndex != i) ? t3.First<Triangle>(p => p.TriangleIndex != i) : null;
        }

        return triangleNeighbors;
    }

    public static Triangle[] FindTriangleNeighbors(Vector3[] vertices, int[] triangles)
    {
        Mesh mesh = new Mesh();
        mesh.vertices = vertices;
        mesh.triangles = triangles;
        return FindTriangleNeighbors(mesh);
    }

    public static Cell[] ComputeBackgroundGrid(Bounds meshBounds, Vector3 cellSize, int nbVertices)
    {
        List<Cell> gridBounds = new List<Cell>();

        int nx = Mathf.CeilToInt(meshBounds.size.x / cellSize.x);
        int ny = Mathf.CeilToInt(meshBounds.size.y / cellSize.y);
        int nz = Mathf.CeilToInt(meshBounds.size.z / cellSize.z);

        const float margin = 0.05f;

        Vector3 newCellSize = new Vector3(((meshBounds.size.x + 2 * margin) / nx), ((meshBounds.size.y + 2 * margin) / ny), ((meshBounds.size.z + 2 * margin) / nz));
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nz; k++)
                {
                    Vector3 center = new Vector3(newCellSize.x * (i + 0.5f) + meshBounds.min.x - margin,
                                                newCellSize.y * (j + 0.5f) + meshBounds.min.y - margin,
                                                newCellSize.z * (k + 0.5f) + meshBounds.min.z - margin);

                    gridBounds.Add(new Cell(new Bounds(center, newCellSize), nbVertices / 2));
                }
            }
        }

        return gridBounds.ToArray();
    }

    //Only tried on sphere (not on concave forms)
    public static (Vector3, float) ComputeVolumeAndBarycentre(Vector3[] vertices, int[] triangles, Transform transform)
    {
        float volume = 0f;
        Vector3 barycentre = Vector3.zero;


        for (int i = 0; i < triangles.Length; i += 3)
        {
            Vector3 tempCrossProd = Vector3.Cross(vertices[triangles[i]] - vertices[triangles[i + 2]], vertices[triangles[i + 1]] - vertices[triangles[i + 2]]);

            float tetraVolume = Mathf.Abs(Vector3.Dot(vertices[triangles[i]], Vector3.Cross(vertices[triangles[i + 1]], vertices[triangles[i + 2]]))) / 6f;
            Vector3 tetraCentroid = (vertices[triangles[i]] + vertices[triangles[i + 1]] + vertices[triangles[i + 2]]) * tetraVolume * 0.25f;

            Vector3 normal = Vector3.Cross(vertices[triangles[i + 1]] - vertices[triangles[i]], vertices[triangles[i + 2]] - vertices[triangles[i]]);
            float faceSign = Mathf.Sign(Vector3.Dot(normal, vertices[triangles[i]]));

            Debug.DrawRay(transform.TransformPoint((vertices[triangles[i]] + vertices[triangles[i + 2]] + vertices[triangles[i + 2]]) / 3f), normal);


            volume += faceSign * tetraVolume;
            barycentre += faceSign * tetraCentroid;

        }

        barycentre /= volume;

        return (barycentre, volume);
    }

    public static (Vector3, float, Vector3, Vector3) ComputeVolumeAndBarycentre(Triangle[] triangles, Vector3 linearSpeed, Vector3 angularSpeed, float C = 1f, float rho = 1000, float mu = 1E-3f)
    {
        float volume = 0f;
        Vector3 barycentre = Vector3.zero;
        Vector3 linearDrag = Vector3.zero;
        Vector3 angularDrag = Vector3.zero;

        for (int i = 0; i < triangles.Length; i++)
        {
            Vector3 vertex1 = triangles[i].Vertex1;
            Vector3 vertex2 = triangles[i].Vertex2;
            Vector3 vertex3 = triangles[i].Vertex3;
            Vector3 trianglePosition = (vertex1 + vertex2 + vertex3) / 3f;

            Vector3 normal = Vector3.Cross(vertex1 - vertex2, vertex1 - vertex3);

            float tetraVolume = Mathf.Abs(Vector3.Dot(vertex1, Vector3.Cross(vertex2, vertex3))) / 6f;
            Vector3 tetraCentroid = (vertex1 + vertex2 + vertex3) * tetraVolume * 0.25f;

            float faceSign = Mathf.Sign(Vector3.Dot(normal, vertex1));

            volume += faceSign * tetraVolume;
            barycentre += faceSign * tetraCentroid;

            Vector3 speedAtTriangle = Vector3.Cross(angularSpeed, trianglePosition);//  / (trianglePosition.magnitude * trianglePosition.magnitude);

            //Friction torque
            float projectedSurfaceAngularSpeed = Vector3.Dot(normal, angularSpeed) > 0 ? ProjectedSurface(vertex1, vertex2, vertex3, speedAtTriangle) : 0;
            angularDrag -= rho * 0.5f * speedAtTriangle.magnitude * speedAtTriangle.magnitude *
                            projectedSurfaceAngularSpeed * trianglePosition.magnitude * angularSpeed.normalized * C;

            //Shear stress
            float projectedSurfaceShearStress = ProjectedSurface(vertex1, vertex2, vertex3, trianglePosition);
            angularDrag -= projectedSurfaceShearStress * mu * angularSpeed * trianglePosition.magnitude;

            //Linear friction drag
            float projectedSurfaceLinearSpeed = Vector3.Dot(normal, linearSpeed) > 0 ? ProjectedSurface(vertex1, vertex2, vertex3, linearSpeed) : 0;
            linearDrag -= rho * 0.5f * linearSpeed.magnitude * linearSpeed * projectedSurfaceLinearSpeed * C;

        }

        barycentre /= volume;
        return (barycentre, volume, linearDrag, angularDrag);
    }

    public static Mesh WeldVertices(Mesh aMesh, float aMaxDelta = 0.01f)
    {
        var verts = aMesh.vertices;
        Dictionary<Vector3, int> duplicateHashTable = new Dictionary<Vector3, int>();
        List<int> newVerts = new List<int>();
        int[] map = new int[verts.Length];

        //create mapping and find duplicates, dictionaries are like hashtables, mean fast
        for (int i = 0; i < verts.Length; i++)
        {
            if (!duplicateHashTable.ContainsKey(verts[i]))
            {
                duplicateHashTable.Add(verts[i], newVerts.Count);
                map[i] = newVerts.Count;
                newVerts.Add(i);
            }
            else
            {
                map[i] = duplicateHashTable[verts[i]];
            }
        }

        // create new vertices
        var verts2 = new Vector3[newVerts.Count];
        var normals2 = new Vector3[newVerts.Count];
        var uvs2 = new Vector2[newVerts.Count];
        for (int i = 0; i < newVerts.Count; i++)
        {
            int a = newVerts[i];
            verts2[i] = verts[a];
        }
        // map the triangle to the new vertices
        var tris = aMesh.triangles;
        for (int i = 0; i < tris.Length; i++)
        {
            tris[i] = map[tris[i]];
        }
        aMesh.triangles = tris;
        aMesh.vertices = verts2;

        aMesh.RecalculateBounds();
        aMesh.RecalculateNormals();

        return aMesh;
    }


    public static (Vector3[], int[], Vector2[]) GenerateGridMesh(float width, int xSize, int ySize)
    {

        Vector3[] vertices = new Vector3[(xSize + 1) * (ySize + 1)];
        Vector2[] uvs = new Vector2[(xSize + 1) * (ySize + 1)];
        for (int i = 0, y = 0; y <= ySize; y++)
        {
            for (int x = 0; x <= xSize; x++, i++)
            {
                vertices[i] = new Vector3(x / (float)xSize * width - width / 2, 0, y / (float)xSize * width - width / 2);
                uvs[i] = new Vector2(x / (float)xSize, y / (float)ySize);
            }
        }

        int[] triangles = new int[xSize * ySize * 6];
        for (int ti = 0, vi = 0, y = 0; y < ySize; y++, vi++)
        {
            for (int x = 0; x < xSize; x++, ti += 6, vi++)
            {
                triangles[ti] = vi;
                triangles[ti + 3] = triangles[ti + 2] = vi + 1;
                triangles[ti + 4] = triangles[ti + 1] = vi + xSize + 1;
                triangles[ti + 5] = vi + xSize + 2;
            }
        }

        return (vertices, triangles, uvs);
    }

    private static float ProjectedSurface(Vector3 v1, Vector3 v2, Vector3 v3, Vector3 normal)
    {
        Vector3 vertex1Projected = v1 - Vector3.Dot(v1, normal.normalized) * normal.normalized;
        Vector3 vertex2Projected = v2 - Vector3.Dot(v2, normal.normalized) * normal.normalized;
        Vector3 vertex3Projected = v3 - Vector3.Dot(v3, normal.normalized) * normal.normalized;
        return (Vector3.Cross(vertex1Projected - vertex2Projected, vertex1Projected - vertex3Projected)).magnitude / 2f;

    }

}
