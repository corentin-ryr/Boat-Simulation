using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class FloatingSphere : MonoBehaviour
{
    Isocahedron isocahedron;
    public ComputeShader computeShader;
    MeshFilter meshFilter;

    MeshCollider meshCollider;

    Floater floater;

    public int nbSubdivision;
    // Start is called before the first frame update
    void Start()
    {
        isocahedron = new Isocahedron(nbSubdivision, 1f, computeShader, Callback);
        isocahedron.CreateSphere();
        meshFilter = GetComponent<MeshFilter>();
        meshCollider = GetComponent<MeshCollider>();
        floater = GetComponent<Floater>();
    }

    void Callback(Vector3Int[] trianglesArray, Vector3[] vertices)
    {
        meshFilter.mesh.vertices = isocahedron.vertices.ToArray();
        List<int> triangles = new List<int>();
        foreach (Vector3Int v3 in trianglesArray)
        {
            triangles.Add(v3.x);
            triangles.Add(v3.y);
            triangles.Add(v3.z);
        }
        meshFilter.mesh.triangles = triangles.ToArray();

        meshCollider.sharedMesh = meshFilter.mesh;

        floater.FloatingMesh = meshFilter.mesh;

    }
}
