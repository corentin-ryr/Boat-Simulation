using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public struct Triangle
{
    public Triangle(int triangleIndex, int n1, int n2, int n3) //Vector3 vertex1, Vector3 vertex2, Vector3 vertex3, 
    {
        this.TriangleIndex = triangleIndex;
        // this.Vertex1 = vertex1;
        // this.Vertex2 = vertex2;
        // this.Vertex3 = vertex3;

        this.N1 = n1;
        this.N2 = n2;
        this.N3 = n3;
    }

    public int TriangleIndex { get; }

    //The 3 vertices of the triangle
    // public Vector3 Vertex1 { get; }
    // public Vector3 Vertex2 { get; }
    // public Vector3 Vertex3 { get; }

    //The 3 neighbors
    public int N1 { get; }//Neighbor between triangle[TriangleIndex * 3] and triangle[TriangleIndex * 3+1]
    public int N2 { get; }//Neighbor between triangle[TriangleIndex * 3+1] and triangle[TriangleIndex * 3+2]
    public int N3 { get; }//Neighbor between triangle[TriangleIndex * 3+2] and triangle[TriangleIndex * 3]

    public int[] GetNeighborIndices()
    {
        return new int[] { N1, N2, N3 };
    }

    public override int GetHashCode()
    {
        return TriangleIndex;
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
    public List<Triangle> TriangleSet1 { get => triangleSet1; }
    bool hasCandidatePotential;
    public bool HasCandidatePotential { get => hasCandidatePotential; }

    public void AddSet1(Triangle triangle)
    {
        triangleSet1.Add(triangle);
        hasCandidatePotential = true;
    }
    private List<Triangle> triangleSet2;
    public List<Triangle> TriangleSet2 { get => triangleSet2; }

    public void AddSet2(Triangle triangle)
    {
        triangleSet2.Add(triangle);
    }

    public bool HasCandidates()
    {
        // Debug.Log(TriangleSet1.Count);
        // Debug.Log(TriangleSet2.Count);
        if (triangleSet1.Count > 0 && triangleSet2.Count > 0)
        {
            return true;
        }
        return false;
    }

    public void resetSet2()
    {
        TriangleSet2.Clear();
    }


}

public static class MeshHelper
{

    public static Triangle[] FindTriangleNeighbors(Mesh mesh)
    {
        int Nt = mesh.triangles.Length / 3;
        int Nv = mesh.vertexCount;

        List<int>[] S = new List<int>[Nv];
        for (int i = 0; i < Nt; i++)
        {
            (S[mesh.triangles[i * 3]] ??= new List<int>()).Add(i);
            (S[mesh.triangles[i * 3 + 1]] ??= new List<int>()).Add(i);
            (S[mesh.triangles[i * 3 + 2]] ??= new List<int>()).Add(i);
        }

        Triangle[] triangleNeighbors = new Triangle[Nt];
        for (int i = 0; i < Nt; i++)
        {
            IEnumerable<int> n1 = S[mesh.triangles[i * 3]].Intersect(S[mesh.triangles[i * 3 + 1]]);
            IEnumerable<int> n2 = S[mesh.triangles[i * 3 + 1]].Intersect(S[mesh.triangles[i * 3 + 2]]);
            IEnumerable<int> n3 = S[mesh.triangles[i * 3 + 2]].Intersect(S[mesh.triangles[i * 3]]);

            triangleNeighbors[i] = new Triangle(i,
                                                n1.Any(p => p != i) ? n1.First<int>(p => p != i) : -1,
                                                n2.Any(p => p != i) ? n2.First<int>(p => p != i) : -1,
                                                n3.Any(p => p != i) ? n3.First<int>(p => p != i) : -1);
        }

        return triangleNeighbors;
    }

    public static Cell[] ComputeBackgroundGrid(Bounds meshBounds, Vector3 cellSize, int nbVertices)
    {
        List<Cell> gridBounds = new List<Cell>();

        int nx = Mathf.CeilToInt(meshBounds.size.x / cellSize.x);
        int ny = Mathf.CeilToInt(meshBounds.size.y / cellSize.y);
        int nz = Mathf.CeilToInt(meshBounds.size.z / cellSize.z);
        Debug.Log(nz);
        Vector3 newCellSize = new Vector3((meshBounds.size.x / nx), (meshBounds.size.y / ny), (meshBounds.size.z / nz));
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nz; k++)
                {
                    Vector3 center = new Vector3(newCellSize.x * i + newCellSize.x / 2 + meshBounds.min.x,
                                                newCellSize.y * j + newCellSize.y / 2 + meshBounds.min.y,
                                                newCellSize.z * k + newCellSize.z / 2 + meshBounds.min.z);

                    gridBounds.Add(new Cell(new Bounds(center, newCellSize), nbVertices / 2));
                }
            }
        }

        return gridBounds.ToArray();
    }

    //Only tried on sphere (not on concave forms)
    public static (Vector3, float) ComputeVolumeAndBarycentre(Vector3[] vertices, int[] triangles)
    {
        float volume = 0f;
        Vector3 barycentre = Vector3.zero;


        for (int i = 0; i < triangles.Length; i += 3)
        {
            Vector3 tempCrossProd = Vector3.Cross(vertices[triangles[i]] - vertices[triangles[i + 2]], vertices[triangles[i + 1]] - vertices[triangles[i + 2]]);

            float tetraVolume = Mathf.Abs(Vector3.Dot(-vertices[triangles[i + 2]], tempCrossProd)) / 6f;
            Vector3 tetraCentroid = (vertices[triangles[i]] + vertices[triangles[i + 1]] + vertices[triangles[i + 2]]) * tetraVolume * 0.25f;

            Vector3 normal = Vector3.Cross(vertices[triangles[i + 1]] - vertices[triangles[i]], vertices[triangles[i + 2]] - vertices[triangles[i]]);
            float faceSign = Mathf.Sign(Vector3.Dot(normal, vertices[triangles[i]]));

            volume += faceSign * tetraVolume;
            barycentre += faceSign * tetraCentroid;

        }

        barycentre /= volume;

        return (barycentre, volume);
    }

}
