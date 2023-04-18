using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class FlatWater : MonoBehaviour, IWater
{

    MeshFilter meshFilter;
    Mesh sharedMeshCollider;
    public Mesh Mesh { get => meshFilter.mesh; }

    float[] heightMap;
    public float[] HeightMap { get => heightMap; }

    [Header("Debug option")]
    public bool showBackgroundGrid;

    Triangle[] triangleNeighbors;
    public Triangle[] TriangleNeighbors { get => triangleNeighbors; }

    [Header("Parameters")]
    public float waterAmplitude;
    public int numGrid = 20;
    public float mapSize = 10;



    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();

        (Vector3[] vertices, int[] triangles, Vector2[] uvs) = MeshHelper.GenerateGridMesh(mapSize, numGrid, numGrid);
        meshFilter.mesh.vertices = vertices;
        meshFilter.mesh.triangles = triangles;
        meshFilter.mesh.uv = uvs;
        meshFilter.mesh.RecalculateNormals();

        heightMap = new float[vertices.Length];

        gameObject.layer = 4; //4 is Water

        // MeshDataPrecomputation();
    }


    void FixedUpdate()
    {
        // Vector3[] vertices = meshFilter.mesh.vertices;

        // for (int i = 0; i < vertices.Length; i++)
        // {
        //     heightMap[i] = GetWaterHeight(vertices[i]);
        //     vertices[i].y = heightMap[i];
        // }

        // for (int i = 0; i < triangleNeighbors.Length; i++)
        // {
        //     triangleNeighbors[i].SetVertex1Height(GetWaterHeight(triangleNeighbors[i].Vertex1));
        //     triangleNeighbors[i].SetVertex2Height(GetWaterHeight(triangleNeighbors[i].Vertex2));
        //     triangleNeighbors[i].SetVertex3Height(GetWaterHeight(triangleNeighbors[i].Vertex3));
        // }

        // meshFilter.mesh.vertices = vertices;
    }

    private void MeshDataPrecomputation()
    {
        (triangleNeighbors, _) = MeshHelper.FindTriangleNeighbors(meshFilter.mesh);
    }

    public float GetWaterHeight(Vector3 position)
    {
        return 0;
    }

    #region Debug and Gizmos =======================================================================
    void DrawTriangleHelper(int i)
    {
        Vector3 v1 = transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[i * 3]]);
        Vector3 v2 = transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[i * 3 + 1]]);
        Vector3 v3 = transform.TransformPoint(meshFilter.mesh.vertices[meshFilter.mesh.triangles[i * 3 + 2]]);

        Debug.DrawRay(Vector3.zero, (v1 + v2 + v3) / 3f, Color.red, 5000f);
    }
    #endregion


}