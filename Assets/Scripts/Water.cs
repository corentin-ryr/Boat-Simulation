using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Water : MonoBehaviour
{

    MeshFilter meshFilter;
    Mesh mesh;
    public Mesh Mesh { get => mesh; }

    [Header("Debug option")]
    public bool showBackgroundGrid;

    Triangle[] triangleNeighbors;
    public Triangle[] TriangleNeighbors { get => triangleNeighbors; }


    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();
        mesh = meshFilter.mesh; ;

        MeshDataPrecomputation();
    }


    void FixedUpdate()
    {
        Vector3[] vertices = mesh.vertices;

        for (int i = 0; i < vertices.Length; i++)
        {
            vertices[i] = new Vector3(vertices[i].x, Mathf.Sin((vertices[i].x + Time.time) * 1f) * 0.2f, vertices[i].z);
        }

        mesh.vertices = vertices;
        meshFilter.mesh = mesh;
    }

    private void MeshDataPrecomputation()
    {
        triangleNeighbors = MeshHelper.FindTriangleNeighbors(mesh);
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
