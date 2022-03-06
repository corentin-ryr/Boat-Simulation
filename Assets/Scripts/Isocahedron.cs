using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;
using UnityEngine.Rendering;


public class Isocahedron
{
    #region Variables
    //Parameters
    private float radius;
    private int nbSubdivision;

    private ComputeShader subdivisionShader;
    private ComputeShader shapeShader;

    //Elements of the unity mesh
    public List<Vector3> vertices = new List<Vector3>();
    public FacesAndEdgesList faces { get; }
    // public Vector3Int[] triangles_array = new Vector3Int[0];
    private NoDuplicatesList edges;

    System.Action<Vector3Int[], Vector3[]> callback;


    #endregion

    public Isocahedron(int nbSubdivision, float radius, ComputeShader subdivisionShader, System.Action<Vector3Int[], Vector3[]> callback )
    {
        this.nbSubdivision = nbSubdivision;
        this.radius = radius;
        this.subdivisionShader = subdivisionShader;

        faces = new FacesAndEdgesList();
        edges = faces.getEdges();
        this.callback = callback;
    }

    public void CreateSphere()
    {
        CreateVertices(); // Creates only the points (vertices)
        TriangulateVertices(); // Create the faces, aka the triangles between the points
        RefineSphereGPU(nbSubdivision); // Subdivide every face by adding new points and triangulating them
    }

    public void SubdivideFace(Vector3[] initialVerticesArray)
    {

        vertices = new List<Vector3>(initialVerticesArray);
        faces.Add(new Vector3Int(2, 1, 0));

        RefineSphereGPU(nbSubdivision);
    }


    private void OnCompleteReadback(AsyncGPUReadbackRequest request, int whichArray, ref ComputeBuffer buffer)
    {
        if (request.hasError)
        {
            UnityEngine.Debug.Log("GPU readback error detected.");
            return;
        }
        if (whichArray == 0)
        {
            vertices = request.GetData<Vector3>().ToList<Vector3>();

        }
        else if (whichArray == 1)
        {
            Vector3Int[] triangles_array = request.GetData<Vector3Int>().ToArray();

            faces.Clear(); //Clears edges too
            foreach (Vector3Int item in triangles_array) //It updates the edges at the same time
            {
                faces.Add(item);
            }

            callback(triangles_array, vertices.ToArray());
        }
        buffer.Dispose();

    }

    #region Mesh generation
    private void CreateVertices()
    {
        // create 12 vertices of a icosahedron
        float t = (float)(1.0 + Math.Sqrt(5.0)) / 2.0f; //Golden ratio

        vertices.Add(new Vector3(-1, t, 0));
        vertices.Add(new Vector3(1, t, 0));
        vertices.Add(new Vector3(-1, -t, 0));
        vertices.Add(new Vector3(1, -t, 0));

        vertices.Add(new Vector3(0, -1, t));
        vertices.Add(new Vector3(0, 1, t));
        vertices.Add(new Vector3(0, -1, -t));
        vertices.Add(new Vector3(0, 1, -t));

        vertices.Add(new Vector3(t, 0, -1));
        vertices.Add(new Vector3(t, 0, 1));
        vertices.Add(new Vector3(-t, 0, -1));
        vertices.Add(new Vector3(-t, 0, 1));
    }

    private void TriangulateVertices()
    {
        // create 20 triangles of the icosahedron

        // 5 faces around point 0
        faces.Add(new Vector3Int(0, 11, 5));
        faces.Add(new Vector3Int(0, 5, 1));
        faces.Add(new Vector3Int(0, 1, 7));
        faces.Add(new Vector3Int(0, 7, 10));
        faces.Add(new Vector3Int(0, 10, 11));

        // 5 adjacent faces
        faces.Add(new Vector3Int(1, 5, 9));
        faces.Add(new Vector3Int(5, 11, 4));
        faces.Add(new Vector3Int(11, 10, 2));
        faces.Add(new Vector3Int(10, 7, 6));
        faces.Add(new Vector3Int(7, 1, 8));

        // 5 faces around point 3
        faces.Add(new Vector3Int(3, 9, 4));
        faces.Add(new Vector3Int(3, 4, 2));
        faces.Add(new Vector3Int(3, 2, 6));
        faces.Add(new Vector3Int(3, 6, 8));
        faces.Add(new Vector3Int(3, 8, 9));

        // 5 adjacent faces
        faces.Add(new Vector3Int(4, 9, 5));
        faces.Add(new Vector3Int(2, 4, 11));
        faces.Add(new Vector3Int(6, 2, 10));
        faces.Add(new Vector3Int(8, 6, 7));
        faces.Add(new Vector3Int(9, 8, 1));
    }


    private void RefineSphereGPU(int nbSubdivision)
    {
        int kernelSubdivide = subdivisionShader.FindKernel("SubdivideEdges");
        int kernelTriangulate = subdivisionShader.FindKernel("CreateFaces");

        // Send the buffers to the Compute Shader ======================================================
        int verticesBufferLength = vertices.Count + (nbSubdivision - 1) * edges.Count + faces.Count * (nbSubdivision - 2) * (nbSubdivision - 1) / 2;
        Vector3[] verticesArray = vertices.ToArray();
        Array.Resize(ref verticesArray, verticesBufferLength);
        ComputeBuffer verticesBuffer = ComputeHelper.CreateAndSetBuffer(verticesArray, subdivisionShader, "vertices", kernelSubdivide);
        subdivisionShader.SetBuffer(kernelTriangulate, "vertices", verticesBuffer);

        ComputeBuffer edgesBuffer = ComputeHelper.CreateAndSetBuffer(edges.ToArray(), subdivisionShader, "edges", kernelSubdivide);

        int lenghtKeys = edges.Count + faces.Count * (nbSubdivision - 1);
        int[] keysArray = new int[lenghtKeys];
        ComputeBuffer keysBuffer = ComputeHelper.CreateAndSetBuffer(keysArray, subdivisionShader, "keys", kernelSubdivide);
        subdivisionShader.SetBuffer(kernelTriangulate, "keys", keysBuffer);

        int[] cacheArray = new int[lenghtKeys * (nbSubdivision + 1)];
        ComputeBuffer cacheBuffer = ComputeHelper.CreateAndSetBuffer(cacheArray, subdivisionShader, "cache", kernelSubdivide);
        subdivisionShader.SetBuffer(kernelTriangulate, "cache", cacheBuffer);

        Vector3Int[] trianglesArray = new Vector3Int[faces.Count * nbSubdivision * nbSubdivision];
        for (int i = 0; i < trianglesArray.Length; i++)
        {
            if (i % (nbSubdivision * nbSubdivision) == 0) { trianglesArray[i] = faces[i / (nbSubdivision * nbSubdivision)]; }
            else { trianglesArray[i] = new Vector3Int(0, 0, 0); }
        }
        ComputeBuffer trianglesBuffer = ComputeHelper.CreateAndSetBuffer(trianglesArray, subdivisionShader, "triangles", kernelTriangulate);

        // Set the int and float values ==================================================================
        subdivisionShader.SetInt("nbVertices", vertices.Count);
        subdivisionShader.SetInt("nbEdges", edges.Count);
        subdivisionShader.SetInt("nbFaces", faces.Count);
        subdivisionShader.SetInt("nbSubdivision", nbSubdivision);
        subdivisionShader.SetFloat("radius", radius);

        // Launch the computation ========================================================================
        ComputeHelper.Run(subdivisionShader, edges.Count, 1, 1, kernelSubdivide);
        ComputeHelper.Run(subdivisionShader, faces.Count, 1, 1, kernelTriangulate);

        // Retrieve the data =============================================================================
        AsyncGPUReadback.Request(verticesBuffer, request => OnCompleteReadback(request, 0, ref verticesBuffer));
        AsyncGPUReadback.Request(trianglesBuffer, request => OnCompleteReadback(request, 1, ref trianglesBuffer));

        // Free the data on the GPU
        edgesBuffer.Dispose();
        keysBuffer.Dispose();
        cacheBuffer.Dispose();
    }
    #endregion


}
