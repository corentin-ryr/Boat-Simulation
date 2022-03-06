using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

[RequireComponent(typeof(MeshCollider))]
[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class Water : MonoBehaviour
{

    MeshFilter meshFilter;
    MeshCollider meshCollider;
    // Mesh mesh;
    Mesh sharedMesh;
    public Mesh Mesh { get => meshFilter.mesh; }

    [Header("Debug option")]
    public bool showBackgroundGrid;

    Triangle[] triangleNeighbors;
    public Triangle[] TriangleNeighbors { get => triangleNeighbors; }


    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();

        (Vector3[] vertices, int[] triangles) = MeshHelper.GenerateGridMesh(20, 20);
        meshFilter.mesh.vertices = vertices;
        meshFilter.mesh.triangles = triangles;

        sharedMesh = new Mesh();
        sharedMesh.vertices = meshFilter.mesh.vertices;
        sharedMesh.triangles = meshFilter.mesh.triangles.Reverse().ToArray();


        meshCollider = GetComponent<MeshCollider>();
        meshCollider.sharedMesh = sharedMesh;


        gameObject.layer = 4; //4 is Water

        MeshDataPrecomputation();
    }


    void FixedUpdate()
    {
        Vector3[] vertices = meshFilter.mesh.vertices;

        for (int i = 0; i < vertices.Length; i++)
        {
            vertices[i].y = GetHeight(vertices[i]);
        }

        for (int i = 0; i < triangleNeighbors.Length; i++)
        {
            triangleNeighbors[i].SetVertex1Height(GetHeight(triangleNeighbors[i].Vertex1));
            triangleNeighbors[i].SetVertex2Height(GetHeight(triangleNeighbors[i].Vertex2));
            triangleNeighbors[i].SetVertex3Height(GetHeight(triangleNeighbors[i].Vertex3));
        }

        // mesh.vertices = vertices;
        meshFilter.mesh.vertices = vertices;

        // sharedMesh.vertices = vertices;
        meshCollider.sharedMesh.vertices = vertices;
    }

    private void MeshDataPrecomputation()
    {
        triangleNeighbors = MeshHelper.FindTriangleNeighbors(meshFilter.mesh);
    }

    private float GetHeight(Vector3 position)
    {
        return Mathf.Sin(position.x * 0.5f + Time.time) * 1f;
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
