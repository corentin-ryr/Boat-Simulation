﻿using UnityEngine;
using System.Collections.Generic;
using System;
using System.Linq;

[RequireComponent(typeof(Rigidbody), typeof(MeshCollider))]
public class Floater : MonoBehaviour
{
    #region Variables =====================================================================================

    //Internal variables of the mesh
    float volume;
    float newVolume;
    Vector3 barycentre;
    Vector3 newBarycentre;
    bool meshReady; //Set after the event

    //Data structures and vectors (global for optimization purposes)
    Triangle[] boatTriangleNeighbors;
    Triangle[] waterTriangleNeighbors;
    Cell[] gridCells;
    Triangle[] triangleCandidateBoat;
    List<Triangle> triangleCandidateWater;

    List<Vector3> chain;
    List<Vector3> lowerChain;
    Vector3 waterIntersectionNormal;
    Vector3 averageWaterPosition;
    int kernelTriangleCandidate;

    ComputeBuffer tempBuffer;
    ComputeBuffer bufferNbTriangleIntersects;
    ComputeBuffer bufferTriangleIntersects;

    MovingAverage smoothedVector;
    Vector3[] waterVertices;
    int[] waterTriangles;

    private float boatMaxSize;

    //References ================
    Rigidbody boatRigidbody;
    WaterManager water;

    //Editor variables ================================
    [Header("References")]
    public ComputeShader triangleCandidateShader;

    [SerializeField]
    private Mesh floatingMesh;

    [Header("Parameters")]
    public Vector3 cellMaxSize = new Vector3(0.2f, 0.2f, 0.2f);
    public Vector3 centerOfMass; // To adjust the center of mass if wanted
    public float C = 1;
    public int waterGridWidth = 5;
    public bool useColliderAsMesh = true; // Wether to get the mesh from the collider and use it as the floating mesh


    [Header("Debug option")]
    public bool showBackgroundGrid;
    public bool showStatValues;
    public bool showCandidatesBoat;
    public bool showCandidatesWater;
    public bool showIntersectionLine;
    public bool showIntersectionMesh;
    public bool showLowerChain;
    public bool showWaterNormal;
    public bool showImmersedBarycentre;
    public bool showBelowWaterPoints;
    public bool showCenterOfMass;
    public bool showWaterMesh;
    public bool showForces;
    public bool showWaterSurface;


    #endregion


    public Mesh FloatingMesh
    {
        get => floatingMesh;
        set
        {
            MeshDataPrecomputation();
            floatingMesh = value;
        }
    }

    // Start is called before the first frame update
    void Start()
    {
        water = GameObject.FindObjectOfType<WaterManager>();
        triangleCandidateShader = Instantiate(triangleCandidateShader);

        boatRigidbody = GetComponent<Rigidbody>();
        boatRigidbody.centerOfMass = centerOfMass;

        if (useColliderAsMesh)
        {
            MeshCollider meshCollider = GetComponent<MeshCollider>();
            floatingMesh = meshCollider.sharedMesh;
        }

        chain = new List<Vector3>();
        lowerChain = new List<Vector3>();
        triangleCandidateWater = new List<Triangle>();

        smoothedVector = new MovingAverage(5);

        if (floatingMesh != null)
        {
            MeshDataPrecomputation();
            meshReady = true;
        }

    }


    void Update()
    {
        if (!meshReady) return;

        if (showStatValues) //Used to check if the volume and barycentre computation is working well
        {
            UnityEngine.Debug.Log("Full volume " + volume);
            UnityEngine.Debug.Log("Geometric center of the boat " + barycentre);
            UnityEngine.Debug.Log("Immersed volume " + newVolume);
            UnityEngine.Debug.Log("Geometric center of the immersed part " + newBarycentre);
        }

        if (showBackgroundGrid)
        {
            UnityEngine.Debug.Log("Number of cells " + gridCells.Length);
            for (int i = 0; i < gridCells.Length; i++)
            {
                DebugHelper.DrawBounds(gridCells[i].Bounds, 0, this.gameObject);
            }
        }
        if (showWaterMesh)
        {
            DebugHelper.ShowMesh(waterVertices, waterTriangles, transform, Color.red);
        }
    }

    void FixedUpdate()
    {
        if (!meshReady) return;
        UpdateWaterMesh();

        int[] cellCandidates = FindTriangleCandidatesGPU();

        if (showCandidatesBoat) DebugHelper.ShowMesh(triangleCandidateBoat, transform, Color.red);
        if (showCandidatesWater) DebugHelper.ShowMesh(triangleCandidateWater.ToArray(), transform, Color.blue);


        //TODO sometimes it bugs
        List<Triangle> triangleRing;
        (chain, triangleRing) = FloaterHelper.ComputeIntersectingLine(gridCells, transform);

        if (chain == null)
        {
            //Test if there is water above
            if (water.GetWaterHeight(transform.position) > transform.position.y)
            {
                boatRigidbody.AddForceAtPosition(volume * 1000 * 9.91f * Vector3.up, transform.TransformPoint(barycentre));
            }
            return;
        }

        //Compute water surface normal (we take 3 points as far away from each other as possible and compute the normal)
        Vector3 a = chain[0];
        Vector3 b = chain[chain.Count / 3];
        Vector3 c = chain[chain.Count / 3 * 2];
        waterIntersectionNormal = Vector3.Cross(a - b, a - c).normalized;
        waterIntersectionNormal = Vector3.Dot(waterIntersectionNormal, transform.InverseTransformPoint(Vector3.up)) < 0 ? -waterIntersectionNormal : waterIntersectionNormal; // We inverse the normal direction if it points downward

        averageWaterPosition = new Vector3(chain.Average(x => x.x), chain.Average(x => x.y), chain.Average(x => x.z));

        //We triangulate the boat (we need to find the triangles below the line)
        Triangle[] bottomHalf = FindBottomHalfOfBoat();
        if (bottomHalf.Length == 0) return;

        //We triangulate the sea (we just create a plane with the circle)
        List<Triangle> seaTriangulated = new List<Triangle>();
        for (int i = 0; i < chain.Count - 1; i++)
        {
            if (Vector3.Dot(Vector3.Cross(averageWaterPosition - chain[i], averageWaterPosition - chain[i + 1]), waterIntersectionNormal) > 0)
            {
                seaTriangulated.Add(new Triangle(0, averageWaterPosition, chain[i], chain[i + 1]));
            }
            else
            {
                seaTriangulated.Add(new Triangle(0, averageWaterPosition, chain[i + 1], chain[i]));
            }
        }

        if (showIntersectionMesh)
        {
            DebugHelper.ShowMesh(bottomHalf, transform, Color.blue, false);
            DebugHelper.ShowMesh(seaTriangulated.ToArray(), transform, Color.blue, false);
        }

        //Apply the forces
        ApplyBuoyancy(bottomHalf, seaTriangulated.ToArray());
    }

    #region Precomputations ======================================================================================
    private void MeshDataPrecomputation() //Executing once the mesh has been created
    {
        floatingMesh = MeshHelper.WeldVertices(floatingMesh, 1E-5f);
        (barycentre, volume) = MeshHelper.ComputeVolumeAndBarycentre(floatingMesh.vertices, floatingMesh.triangles, transform);

        //Precomputing the triangle neighbor relations and the background grid (only once at the beginning because we suppose that the mesh does not change)
        boatTriangleNeighbors = MeshHelper.FindTriangleNeighbors(floatingMesh);

        gridCells = MeshHelper.ComputeBackgroundGrid(floatingMesh.bounds, cellMaxSize, floatingMesh.triangles.Length);

        AssignTrianglesToBoundingBox();

        CreateWaterGrid();
        waterTriangleNeighbors = MeshHelper.FindTriangleNeighbors(waterVertices, waterTriangles);

        //Prepare the compute shader =======================================
        kernelTriangleCandidate = triangleCandidateShader.FindKernel("FindTriangleCandidate");

        //Set the bounds of gridcells
        Vector3[] cellPositions = new Vector3[gridCells.Length];
        Vector3[] cellSizes = new Vector3[gridCells.Length];
        for (int i = 0; i < gridCells.Length; i++)
        {
            cellPositions[i] = gridCells[i].Bounds.center;
            cellSizes[i] = gridCells[i].Bounds.size;
        }
        List<ComputeBuffer> tempBuffers = new List<ComputeBuffer>(); // We don't release those buffers

        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(cellPositions, triangleCandidateShader, "cellPositions", kernelTriangleCandidate));
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(cellSizes, triangleCandidateShader, "cellSizes", kernelTriangleCandidate));
        triangleCandidateShader.SetInt("nbCells", gridCells.Length);

        //Set the water mesh data
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(waterVertices, triangleCandidateShader, "vertices", kernelTriangleCandidate));
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(waterTriangles, triangleCandidateShader, "triangles", kernelTriangleCandidate));
        triangleCandidateShader.SetInt("nbTriangles", waterTriangles.Length / 3);

        meshReady = true;
    }

    private void CreateWaterGrid()
    {
        Bounds bound = floatingMesh.bounds;
        boatMaxSize = (bound.max - bound.center).magnitude * 2;
        (waterVertices, waterTriangles, _) = MeshHelper.GenerateGridMesh(boatMaxSize, 5, 5);
    }

    private void AssignTrianglesToBoundingBox()
    {
        Vector3 a, b, c;
        Vector3 mins, maxes;
        Bounds triangleBoundingbox;
        for (int i = 0; i < floatingMesh.triangles.Length; i += 3)
        {
            a = floatingMesh.vertices[floatingMesh.triangles[i]];
            b = floatingMesh.vertices[floatingMesh.triangles[i + 1]];
            c = floatingMesh.vertices[floatingMesh.triangles[i + 2]];

            mins = Vector3.Min(Vector3.Min(a, b), c);
            maxes = Vector3.Max(Vector3.Max(a, b), c);

            triangleBoundingbox = new Bounds((maxes + mins) / 2, maxes - mins);

            for (int j = 0; j < gridCells.Length; j++)
            {
                if (gridCells[j].Intersects(triangleBoundingbox))
                {
                    gridCells[j].AddSet1(boatTriangleNeighbors[i / 3]);
                }
            }
        }

        //Remove cells with no potential
        List<Cell> newCells = new List<Cell>();
        foreach (Cell cell in gridCells)
        {
            if (cell.HasCandidatePotential) newCells.Add(cell);
        }
        gridCells = newCells.ToArray();
    }
    #endregion

    #region Intersection computation with water ==================================================================

    private int[] FindTriangleCandidates() //Return the indices of the cells with candidate triangle (for finding the first intersection of the chain)
    {
        //Assign water triangles to boat bounding boxes 
        List<int> candidateCells = new List<int>();

        foreach (Cell cell in gridCells)
        {
            cell.resetSet2();
        }

        Vector3 mins, maxes;
        Bounds triangleBoundingbox = new Bounds();
        for (int i = 0; i < waterTriangles.Length; i += 3)
        {
            Vector3 a = waterVertices[waterTriangles[i]];
            Vector3 b = waterVertices[waterTriangles[i + 1]];
            Vector3 c = waterVertices[waterTriangles[i + 2]];

            mins = Vector3.Min(Vector3.Min(a, b), c);
            maxes = Vector3.Max(Vector3.Max(a, b), c);

            triangleBoundingbox.center = (maxes + mins) / 2;//new Bounds((maxes + mins) / 2, maxes - mins);
            triangleBoundingbox.size = maxes - mins;

            for (int j = 0; j < gridCells.Length; j++)
            {
                if (gridCells[j].HasCandidatePotential && gridCells[j].Intersects(triangleBoundingbox))
                {
                    gridCells[j].AddSet2(waterTriangleNeighbors[i / 3]);
                    candidateCells.Add(j);
                }
            }
        }

        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>(floatingMesh.triangles.Length / 3);
        if (showCandidatesWater) triangleCandidateWater.Clear();

        foreach (int cellIndex in candidateCells)
        {
            triangleCandidateBoatHS.UnionWith(gridCells[cellIndex].TriangleSet1);
            if (showCandidatesWater) triangleCandidateWater.AddRange(gridCells[cellIndex].TriangleSet2);
        }

        triangleCandidateBoat = triangleCandidateBoatHS.ToArray();
        return candidateCells.ToArray();
    }

    private int[] FindTriangleCandidatesGPU() //Optimized version of the cpu one above
    {
        //Set the frame data
        int triangleCount = waterTriangles.Length / 3;
        tempBuffer = ComputeHelper.CreateAndSetBuffer(waterVertices, triangleCandidateShader, "vertices", kernelTriangleCandidate);
        bufferNbTriangleIntersects = ComputeHelper.CreateAndSetBuffer<int>(gridCells.Length, triangleCandidateShader, "nbTriangleIntersects", kernelTriangleCandidate);
        bufferTriangleIntersects = ComputeHelper.CreateAndSetBuffer<int>(triangleCount * gridCells.Length, triangleCandidateShader, "triangleIntersects", kernelTriangleCandidate);

        ComputeHelper.Run(triangleCandidateShader, gridCells.Length, 1, 1, kernelTriangleCandidate);

        int[] triangleIntersects = new int[triangleCount * gridCells.Length];
        int[] nbTriangleIntersects = new int[gridCells.Length];

        bufferNbTriangleIntersects.GetData(nbTriangleIntersects); // Freezes when after the other get data (the last get data frozes ?)
        bufferTriangleIntersects.GetData(triangleIntersects); // Freezing line

        bufferTriangleIntersects.Release();
        bufferNbTriangleIntersects.Release();
        tempBuffer.Release();


        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>(floatingMesh.triangles.Length / 3);
        if (showCandidatesWater) triangleCandidateWater.Clear();

        List<int> cellsWithPotential = new List<int>();
        for (int i = 0; i < gridCells.Length; i++)
        {
            gridCells[i].resetSet2();
            if (nbTriangleIntersects[i] > 0)
            {
                cellsWithPotential.Add(i);
                triangleCandidateBoatHS.UnionWith(gridCells[i].TriangleSet1);
                for (int j = 0; j < nbTriangleIntersects[i]; j++)
                {
                    Triangle triangle = waterTriangleNeighbors[triangleIntersects[i * triangleCount + j]];
                    gridCells[i].AddSet2(triangle);
                    if (showCandidatesWater) triangleCandidateWater.Add(triangle);
                }
            }
        }

        triangleCandidateBoat = triangleCandidateBoatHS.ToArray();
        return cellsWithPotential.ToArray();
    }

    private void UpdateWaterMesh()
    {
        for (int i = 0, y = 0; y <= waterGridWidth; y++)
        {
            for (int x = 0; x <= waterGridWidth; x++, i++)
            {
                waterVertices[i] = new Vector3(x / (float)waterGridWidth * boatMaxSize - boatMaxSize / 2, 0, y / (float)waterGridWidth * boatMaxSize - boatMaxSize / 2) + transform.position; //World coordinates
                waterVertices[i] = transform.InverseTransformPoint(new Vector3(waterVertices[i].x, water.GetWaterHeight(waterVertices[i]), waterVertices[i].z));
            }
        }

        for (int i = 0; i < waterTriangleNeighbors.Length; i++)
        {
            waterTriangleNeighbors[i].Vertex1 = waterVertices[waterTriangles[waterTriangleNeighbors[i].N1]];
            waterTriangleNeighbors[i].Vertex2 = waterVertices[waterTriangles[waterTriangleNeighbors[i].N2]];
            waterTriangleNeighbors[i].Vertex3 = waterVertices[waterTriangles[waterTriangleNeighbors[i].N3]];
        }
    }
    #endregion


    #region Triangulate the mesh =================================================================================

    private Triangle[] FindBottomHalfOfBoat()
    {
        float waterSurfaceDepth = PointHeight(averageWaterPosition);

        HashSet<Triangle> bottomHalf = new HashSet<Triangle>();
        lowerChain.Clear();

        //Find a point underwater
        Triangle currentTriangle = triangleCandidateBoat[0];
        float minDepth = float.MaxValue;
        while (minDepth > waterSurfaceDepth)
        {
            minDepth = MaxTriangleHeight(currentTriangle);
            float previousDepth = minDepth;
            //Get the neighbor further down
            foreach (Triangle neighbor in currentTriangle.GetNeighbors())
            {
                float depth = MaxTriangleHeight(neighbor);
                if (depth < minDepth)
                {
                    currentTriangle = neighbor;
                    minDepth = depth;
                }
            }

            if (previousDepth == minDepth) break;
        }

        if (MaxTriangleHeight(currentTriangle) > waterSurfaceDepth) throw new Exception("No triangle underwater");


        //We have an initial triangle underwater
        Queue<Triangle> Q = new Queue<Triangle>();
        Q.Enqueue(currentTriangle);
        bottomHalf.Add(currentTriangle);

        while (Q.Count > 0)
        {
            currentTriangle = Q.Dequeue();
            foreach (Triangle neighbor in currentTriangle.GetNeighbors())
            {
                if (MaxTriangleHeight(neighbor) < waterSurfaceDepth)
                {
                    if (bottomHalf.Add(neighbor)) Q.Enqueue(neighbor);
                }
                else
                {
                    lowerChain.AddRange(currentTriangle.GetVertices().Intersect(neighbor.GetVertices()));
                }
            }
        }

        return bottomHalf.ToArray();
    }

    private Triangle[] TriangulateIntermediate()
    {
        throw new NotImplementedException();
    }

    #endregion

    #region Compute forces on rigidbody ==============================================================================

    private void ApplyBuoyancy(Triangle[] bottomHalf, Triangle[] trianglesSea)
    {
        Vector3 localVelocity = transform.InverseTransformDirection(boatRigidbody.velocity);
        Vector3 localAngularVelocity = transform.InverseTransformDirection(boatRigidbody.angularVelocity);

        (Vector3 barycentreBoat, float volumeBoat, Vector3 linearDragBoat, Vector3 angularDragBoat) = MeshHelper.ComputeVolumeAndBarycentre(bottomHalf, localVelocity, localAngularVelocity, C);
        (Vector3 barycentreSea, float volumeSea, Vector3 linearDragSea, Vector3 angularDragSea) = MeshHelper.ComputeVolumeAndBarycentre(trianglesSea, localVelocity, localAngularVelocity, C, 1.2f, 1E-5f);

        newVolume = volumeSea + volumeBoat;
        newBarycentre = (barycentreBoat * volumeBoat + barycentreSea * volumeSea) / newVolume;

        boatRigidbody.AddForceAtPosition(newVolume * 1000 * 9.81f * transform.TransformDirection(waterIntersectionNormal), transform.TransformPoint(newBarycentre));

        boatRigidbody.AddForce(transform.TransformDirection(linearDragBoat + linearDragSea));
        boatRigidbody.AddTorque(transform.TransformDirection(angularDragSea + angularDragBoat));

        if (newVolume > 4) UnityEngine.Debug.LogWarning("Warning volume" + newVolume);
        if (newBarycentre.magnitude > 1.2) UnityEngine.Debug.LogWarning("Warning barycenter" + newBarycentre);
        if (showWaterNormal) UnityEngine.Debug.DrawRay(transform.position, transform.TransformDirection(waterIntersectionNormal), Color.red);
        if (showForces)
        {
            UnityEngine.Debug.DrawRay(transform.TransformPoint(newBarycentre), newVolume * 9.81f * transform.TransformDirection(waterIntersectionNormal) / 10, Color.green);
            UnityEngine.Debug.DrawRay(transform.position, transform.TransformDirection(linearDragBoat + linearDragSea) / 10, Color.cyan);
        }
    }

    #endregion


    #region Helper functions ===================================================================================
    private float MaxTriangleHeight(Triangle triangle)
    {
        return Mathf.Max(Mathf.Max(
                        PointHeight(triangle.Vertex1),
                        PointHeight(triangle.Vertex2)),
                        PointHeight(triangle.Vertex3));
    }


    /// <summary>
    ///     Compute the height of the point projected on the water normal 
    /// </summary>
    /// <param name="point">In world coordinates</param>
    /// <returns></returns>
    private float PointHeight(Vector3 point)
    {
        return Vector3.Dot(point, waterIntersectionNormal);
    }
    #endregion



    #region Gizmos and debug =====================================================================================

    void OnDrawGizmos()
    {
        if (showStatValues)
        {
            Gizmos.color = Color.red;
            // Gizmos.DrawSphere(transform.TransformPoint(barycentre), 0.1f);
        }

        if (chain?.Count > 0 && showIntersectionLine)
        {
            Gizmos.color = Color.cyan;
            Gizmos.DrawSphere(transform.TransformPoint(chain[0]), 0.01f);
            for (int i = 0; i < chain.Count - 1; i++)
            {
                Gizmos.DrawSphere(transform.TransformPoint(chain[i + 1]), 0.01f);
                Gizmos.DrawLine(transform.TransformPoint(chain[i]), transform.TransformPoint(chain[i + 1]));
            }

        }
        if (lowerChain?.Count > 0 && showLowerChain)
        {
            Gizmos.color = Color.black;
            Gizmos.DrawSphere(transform.TransformPoint(lowerChain[0]), 0.01f);
            for (int i = 0; i < lowerChain.Count - 1; i++)
            {
                Gizmos.DrawSphere(transform.TransformPoint(lowerChain[i + 1]), 0.05f);
            }

        }

        if (showImmersedBarycentre)
        {
            Gizmos.color = Color.red;
            Gizmos.DrawSphere(transform.TransformPoint(newBarycentre), 0.1f);
        }

        if (showBelowWaterPoints)
        {
            Gizmos.color = Color.red;
            Gizmos.DrawSphere(transform.TransformPoint(averageWaterPosition), 0.1f);
            foreach (Vector3 vertex in floatingMesh.vertices)
            {
                if (PointHeight(vertex) < PointHeight(averageWaterPosition))
                {
                    Gizmos.color = Color.blue;
                    Gizmos.DrawSphere(transform.TransformPoint(vertex), 0.01f);
                }
            }
        }

        if (showCenterOfMass && boatRigidbody)
        {
            Gizmos.color = Color.green;
            Gizmos.DrawSphere(transform.TransformPoint(boatRigidbody.centerOfMass), 0.1f);
        }

        if (showWaterSurface)
        {
            DebugHelper.ShowMesh(waterVertices, waterTriangles, gameObject.transform, Color.cyan);
        
        }

    }

    #endregion
}