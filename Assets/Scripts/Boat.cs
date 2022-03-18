using UnityEngine;
using UnityEditor;
using System.Collections.Generic;
using System;
using System.Collections;
using System.Diagnostics;
using System.Linq;
using Habrador_Computational_Geometry;

[RequireComponent(typeof(Rigidbody), typeof(MeshCollider))]
public class Boat : MonoBehaviour
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
    Vector3 a, b, c, d, e, f;
    Vector3 waterIntersectionNormal;
    Vector3 averageWaterPosition;
    int kernelTriangleCandidate;

    ComputeBuffer tempBuffer;
    ComputeBuffer bufferNbTriangleIntersects;
    ComputeBuffer bufferTriangleIntersects;

    MovingAverage smoothedVector;

    //References ================
    MeshFilter meshFilter;
    Rigidbody rigidbody;
    MeshCollider meshCollider;
    Stopwatch stopWatch;

    //Editor variables ============
    [Header("References")]
    public Water water;
    public ComputeShader triangleCandidateShader;

    [Header("Parameters")]
    public Vector3 cellMaxSize = new Vector3(0.2f, 0.2f, 0.2f);
    public bool procedurallyGeneratedMesh;
    public Vector3 centerOfMass;
    public float C = 1;


    [Header("Debug option")]
    public bool showBackgroundGrid;
    public bool showStatValues;
    public bool showCandidatesBoat;
    public bool showCandidatesWater;
    public bool showIntersectionLine;
    public bool showIntersectionMesh;
    public bool showLowerChain;
    public bool showWaterNormal;
    public bool showWaterAboveCheck;
    public bool showImmersedBarycentre;
    public bool showBelowWaterPoints;
    public bool showCenterOfMass;




    #endregion

    void OnEnable()
    {
        EventManager.FinishedGeneratingMesh += MeshDataPrecomputation;
    }
    void OnDisable()
    {
        EventManager.FinishedGeneratingMesh -= MeshDataPrecomputation;
    }

    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();

        rigidbody = GetComponent<Rigidbody>();
        rigidbody.centerOfMass = centerOfMass;

        meshCollider = GetComponent<MeshCollider>();
        stopWatch = new Stopwatch();

        chain = new List<Vector3>();
        lowerChain = new List<Vector3>();
        triangleCandidateWater = new List<Triangle>();

        smoothedVector = new MovingAverage(5);

        Physics.IgnoreCollision(water.GetComponent<Collider>(), meshCollider);

        a = new Vector3();
        b = new Vector3();
        c = new Vector3();
        d = new Vector3();
        e = new Vector3();
        f = new Vector3();

        if (!procedurallyGeneratedMesh)
        {
            EventManager.OnFinishedGeneratingMesh();
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
    }

    void FixedUpdate()
    {
        if (!meshReady) return;
        waterTriangleNeighbors ??= water.TriangleNeighbors;

        int[] cellCandidates = FindTriangleCandidatesGPU();

        if (showCandidatesBoat)
        {
            DebugHelper.ShowMesh(triangleCandidateBoat, transform, Color.red);
        }
        if (showCandidatesWater)
        {
            DebugHelper.ShowMesh(triangleCandidateWater.ToArray(), water.transform, Color.blue);
        }

        //TODO sometimes it bugs
        bool hasIntersection = ComputeIntersectingLine(cellCandidates);

        if (!hasIntersection)
        {
            //Test if there is water above
            if (showWaterAboveCheck) UnityEngine.Debug.DrawRay(transform.position, Vector3.up * 10);
            if (Physics.Raycast(transform.position, Vector3.up, Mathf.Infinity, LayerMask.GetMask("Water")))
            {
                rigidbody.AddForceAtPosition(volume * 1000 * 9.91f * Vector3.up, transform.TransformPoint(barycentre));
            }

            return;
        }

        //Compute water surface normal (we take 3 points as far away from each other as possible and compute the normal)
        a = chain[0];
        b = chain[chain.Count / 3];
        c = chain[chain.Count / 3 * 2];
        waterIntersectionNormal = Vector3.Cross(a - b, a - c).normalized;
        waterIntersectionNormal = Vector3.Dot(waterIntersectionNormal, Vector3.up) < 0 ? -waterIntersectionNormal : waterIntersectionNormal;

        averageWaterPosition = Vector3.zero;
        foreach (Vector3 vector in chain)
        {
            averageWaterPosition += vector;
        }
        averageWaterPosition /= chain.Count;

        if (showWaterNormal)
        {
            UnityEngine.Debug.DrawRay(transform.position, waterIntersectionNormal, Color.red);
        }


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
            DebugHelper.ShowMesh(bottomHalf, transform, Color.blue);
            // DebugHelper.ShowMesh(seaTriangulated.ToArray(), Color.blue);
        }

        //Apply the forces
        ApplyBuoyancy(bottomHalf, seaTriangulated.ToArray());
    }

    #region Precomputations ======================================================================================
    private void MeshDataPrecomputation() //Executing once the mesh has been created
    {
        meshFilter.mesh = MeshHelper.WeldVertices(meshFilter.mesh, 1E-5f);
        (barycentre, volume) = MeshHelper.ComputeVolumeAndBarycentre(meshFilter.mesh.vertices, meshFilter.mesh.triangles, transform);

        meshCollider.sharedMesh = meshFilter.mesh;

        //Precomputing the triangle neighbor relations and the background grid (only once at the beginning because we suppose that the mesh does not change)
        boatTriangleNeighbors = MeshHelper.FindTriangleNeighbors(meshFilter.mesh);

        gridCells = MeshHelper.ComputeBackgroundGrid(meshFilter.mesh.bounds, cellMaxSize, meshFilter.mesh.triangles.Length);

        AssignTrianglesToBoundingBox();

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
        List<ComputeBuffer> tempBuffers = new List<ComputeBuffer>();//We don't release those buffers

        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(cellPositions, triangleCandidateShader, "cellPositions", kernelTriangleCandidate));
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(cellSizes, triangleCandidateShader, "cellSizes", kernelTriangleCandidate));
        triangleCandidateShader.SetInt("nbCells", gridCells.Length);

        //Set the water mesh data
        float[] heightMap = new float[water.Mesh.vertexCount];
        Vector2[] waterVerticesPlanPosition = new Vector2[water.Mesh.vertexCount];
        for (int i = 0; i < water.Mesh.vertexCount; i++)
        {
            heightMap[i] = water.Mesh.vertices[i].y;
            waterVerticesPlanPosition[i] = new Vector2(water.Mesh.vertices[i].x, water.Mesh.vertices[i].z);
        }
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(heightMap, triangleCandidateShader, "verticesHeightMap", kernelTriangleCandidate));
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(waterVerticesPlanPosition, triangleCandidateShader, "verticesPlanPosition", kernelTriangleCandidate));
        tempBuffers.Add(ComputeHelper.CreateAndSetBuffer(water.Mesh.triangles, triangleCandidateShader, "triangles", kernelTriangleCandidate));
        triangleCandidateShader.SetInt("nbTriangles", water.Mesh.triangles.Length / 3);

        //Reset the rigidbody
        rigidbody.inertiaTensor = new Vector3(874, 874, 874);

        meshReady = true;
    }

    private void AssignTrianglesToBoundingBox()
    {
        Vector3 a, b, c;
        Vector3 mins, maxes;
        Bounds triangleBoundingbox;
        for (int i = 0; i < meshFilter.mesh.triangles.Length; i += 3)
        {
            a = meshFilter.mesh.vertices[meshFilter.mesh.triangles[i]];
            b = meshFilter.mesh.vertices[meshFilter.mesh.triangles[i + 1]];
            c = meshFilter.mesh.vertices[meshFilter.mesh.triangles[i + 2]];

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
            if (cell.HasCandidatePotential)
            {
                newCells.Add(cell);
            }
        }
        gridCells = newCells.ToArray();
    }
    #endregion

    #region Intersection computation with water ======================================================

    private int[] FindTriangleCandidates() //Return the indices of the cells with candidate triangle (for finding the first intersection of the chain)
    {
        //Assign water triangles to boat bounding boxes
        Vector3[] vertices = water.Mesh.vertices;
        int[] triangles = water.Mesh.triangles;
        List<int> candidateCells = new List<int>();

        foreach (Cell cell in gridCells)
        {
            cell.resetSet2();
        }

        Matrix4x4 transformation = transform.worldToLocalMatrix * water.transform.localToWorldMatrix;
        Vector3 mins, maxes;

        Bounds triangleBoundingbox = new Bounds();
        for (int i = 0; i < triangles.Length; i += 3)
        {
            a = transformation.MultiplyPoint3x4(vertices[triangles[i]]);
            b = transformation.MultiplyPoint3x4(vertices[triangles[i + 1]]);
            c = transformation.MultiplyPoint3x4(vertices[triangles[i + 2]]);

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

        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>(meshFilter.mesh.triangles.Length / 3);
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
        int triangleCount = water.Mesh.triangles.Length / 3;
        tempBuffer = ComputeHelper.CreateAndSetBuffer(water.HeightMap, triangleCandidateShader, "verticesHeightMap", kernelTriangleCandidate);
        bufferNbTriangleIntersects = ComputeHelper.CreateAndSetBuffer<int>(gridCells.Length, triangleCandidateShader, "nbTriangleIntersects", kernelTriangleCandidate);
        bufferTriangleIntersects = ComputeHelper.CreateAndSetBuffer<int>(triangleCount * gridCells.Length, triangleCandidateShader, "triangleIntersects", kernelTriangleCandidate);

        Matrix4x4 transformation = transform.worldToLocalMatrix * water.transform.localToWorldMatrix;
        triangleCandidateShader.SetMatrix("transformationMatrix", transformation);

        ComputeHelper.Run(triangleCandidateShader, gridCells.Length, 1, 1, kernelTriangleCandidate);

        int[] nbTriangleIntersects = new int[gridCells.Length];
        bufferNbTriangleIntersects.GetData(nbTriangleIntersects);
        bufferNbTriangleIntersects.Release();

        int[] triangleIntersects = new int[gridCells.Length * triangleCount];
        bufferTriangleIntersects.GetData(triangleIntersects);
        bufferTriangleIntersects.Release();

        tempBuffer.Release();

        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>(meshFilter.mesh.triangles.Length / 3);
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

    //The triangle of the boat mesh and the new polygon that will be later triangulated
    Dictionary<Triangle, List<Vector3>> polygons = new Dictionary<Triangle, List<Vector3>>();

    private bool ComputeIntersectingLine(int[] cellCandidates)
    {
        chain.Clear();

        //Find first intersection
        bool swaped = false;
        Triangle currentTriangle;
        Triangle triangleToCheckAgainst = null;
        Vector3 currentP = Vector3.positiveInfinity;
        Triangle[] nextCurrentTriangles = new Triangle[0];
        int currentCase;

        foreach (int cellIndex in cellCandidates)
        {
            Cell cell = gridCells[cellIndex];
            foreach (Triangle triangleSet2 in cell.TriangleSet2)
            {
                foreach (Triangle triangleSet1 in cell.TriangleSet1)
                {
                    (currentCase, nextCurrentTriangles) = IdentifyCase(triangleSet1, meshFilter.mesh, triangleSet2, water.Mesh, ref currentP, ref swaped, ref triangleToCheckAgainst);
                    if (currentCase > 0)
                    {
                        chain.Add(currentP);

                        //Polygon creation
                        if (currentCase == 1)
                        {
                            HandlePolygons(triangleSet1, triangleSet2, currentP, PolygonState.middle);
                        }
                        else
                        {
                            HandlePolygons(nextCurrentTriangles[0], triangleSet2, currentP, PolygonState.start);
                        }
                        goto End;
                    }
                }
            }
        }

        return false;
    End:

        //We have the first triangle, we start looping
        for (int i = 0; i < 50; i++) //Max number of iteration (security)
        {
            foreach (Triangle nextCurrentTriangle in nextCurrentTriangles) //Normally only one triangle in this array but if we have cases 3 or 4 we can have 0 or more than 1 triangles.
            {
                currentTriangle = nextCurrentTriangle;
                (currentCase, nextCurrentTriangles) = IdentifyCase(currentTriangle, swaped ? water.Mesh : meshFilter.mesh, triangleToCheckAgainst,
                                                                    swaped ? meshFilter.mesh : water.Mesh, ref currentP, ref swaped, ref triangleToCheckAgainst);
                if (currentCase > 0)
                {
                    chain.Add(currentP);
                    //Polygon creation
                    if (currentCase == 1)
                    {
                        HandlePolygons(swaped ? triangleToCheckAgainst : currentTriangle, swaped ? currentTriangle : triangleToCheckAgainst, currentP, PolygonState.middle);
                    }
                    else
                    {
                        HandlePolygons(swaped ? triangleToCheckAgainst : currentTriangle, swaped ? currentTriangle : triangleToCheckAgainst, currentP, PolygonState.end);
                        HandlePolygons(nextCurrentTriangles[0], swaped ? currentTriangle : triangleToCheckAgainst, currentP, PolygonState.start);
                    }

                    break;
                }
                UnityEngine.Debug.Log("No intersection");
                return false;
            }

            if ((currentP - chain[0]).magnitude < 1E-3) break;
        }

        return true;
    }

    //Triangle1 is from the boat and triangle2 is from the water
    //TODO compute cases 3 and 4
    private (int, Triangle[]) IdentifyCase(Triangle triangle1, Mesh mesh1, Triangle triangle2, Mesh mesh2, ref Vector3 P, ref bool swaped, ref Triangle triangleToCheckAgainst)
    {
        Matrix4x4 transformation = mesh1 == meshFilter.mesh ? transform.localToWorldMatrix : water.transform.localToWorldMatrix;
        a = transformation.MultiplyPoint3x4(triangle1.Vertex1);
        b = transformation.MultiplyPoint3x4(triangle1.Vertex2);
        c = transformation.MultiplyPoint3x4(triangle1.Vertex3);

        transformation = mesh2 == meshFilter.mesh ? transform.localToWorldMatrix : water.transform.localToWorldMatrix;
        d = transformation.MultiplyPoint3x4(triangle2.Vertex1);
        e = transformation.MultiplyPoint3x4(triangle2.Vertex2);
        f = transformation.MultiplyPoint3x4(triangle2.Vertex3);

        //Case 1 if p is inside of triangle1 
        Triangle[] triangleIndicesToCheck1 = null;
        Vector3 tempP1 = Vector3.negativeInfinity;
        bool case1 = false;
        if (IntersectionEdgeTriangle(a, b, c, d, e, ref tempP1, P))
        {
            case1 = true;
            triangleIndicesToCheck1 = new Triangle[] { triangle2.T1 };
        }
        else if (IntersectionEdgeTriangle(a, b, c, e, f, ref tempP1, P))
        {
            case1 = true;
            triangleIndicesToCheck1 = new Triangle[] { triangle2.T2 };
        }
        else if (IntersectionEdgeTriangle(a, b, c, f, d, ref tempP1, P))
        {
            case1 = true;
            triangleIndicesToCheck1 = new Triangle[] { triangle2.T3 };
        }


        //Case 2 if p is inside of triangle2
        Triangle[] triangleIndicesToCheck2 = null;
        Vector3 tempP2 = Vector3.negativeInfinity;
        bool case2 = false;
        if (IntersectionEdgeTriangle(d, e, f, a, b, ref tempP2, P))
        {
            case2 = true;
            triangleIndicesToCheck2 = new Triangle[] { triangle1.T1 };
        }
        else if (IntersectionEdgeTriangle(d, e, f, b, c, ref tempP2, P))
        {
            case2 = true;
            triangleIndicesToCheck2 = new Triangle[] { triangle1.T2 };
        }
        else if (IntersectionEdgeTriangle(d, e, f, c, a, ref tempP2, P))
        {
            case2 = true;
            triangleIndicesToCheck2 = new Triangle[] { triangle1.T3 };
        }


        if (case1 && case2 && tempP1 == tempP2) //Case when we have the point on the edge of the two triangles. The next triangle to check 
        {
            // triangleToCheck = swaped ? triangleIndicesToCheck2[0] : triangleIndicesToCheck1[0];
            triangleToCheckAgainst = triangleIndicesToCheck2[0];
            swaped = !swaped;
            P = tempP1;
            return (1, triangleIndicesToCheck1); //We act like a case 1 but we have a special triangleToCheck
        }
        else if (case1) //Case 1, we have to swap and we set the triangle to check next
        {
            swaped = !swaped;
            triangleToCheckAgainst = triangle1;
            P = tempP1;
            return (1, triangleIndicesToCheck1);
        }
        else if (case2)//Case 2, we don't swap and we also set the next triangle to check
        {
            P = tempP2;
            triangleToCheckAgainst = triangle2;
            return (2, triangleIndicesToCheck2);
        }

        return (0, triangleIndicesToCheck1);
    }

    private bool IntersectionEdgeTriangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 e, ref Vector3 newP, Vector3 previousP)
    {
        Vector3 n = Vector3.Cross(a - b, a - c);

        float E = Vector3.Dot(a - e, n);
        float D = Vector3.Dot(a - d, n);
        //Are d and e on different side of triangle abc ?
        if (D * E > 0) return false;

        float t = D / (D - E);
        Vector3 p = t * e + (1 - t) * d;

        if (Vector3.Dot(Vector3.Cross(a - b, a - p), n) < -1E-9 ||
            Vector3.Dot(Vector3.Cross(b - c, b - p), n) < -1E-9 ||
            Vector3.Dot(Vector3.Cross(c - a, c - p), n) < -1E-9) return false;

        if (p == previousP)
        {
            return false;
        }
        newP = p;
        return true;
    }

    enum PolygonState
    {
        start,
        middle,
        end
    }
    private void HandlePolygons(Triangle triangle, Triangle waterTriangle, Vector3 newPoint, PolygonState state)
    {
        if (polygons.ContainsKey(triangle))
        {
            if (state == PolygonState.start)
            {
                //Get the point underwater on the intersected edge
                polygons[triangle].Add(PointAdjacentOnEdge(triangle, waterTriangle, newPoint));
            }
            polygons[triangle].Add(newPoint);
            if (state == PolygonState.end)
            {
                //Get the point underwater on the intersected edge
                polygons[triangle].Add(PointAdjacentOnEdge(triangle, waterTriangle, newPoint));
            }
        }
        else if (state == PolygonState.start)
        {
            List<Vector3> tempList = new List<Vector3>();
            tempList.Add(PointAdjacentOnEdge(triangle, waterTriangle, newPoint));
            tempList.Add(newPoint);
            polygons.Add(triangle, tempList);
        }
    }

    private Vector3 PointAdjacentOnEdge(Triangle triangle, Triangle waterTriangle, Vector3 point)
    {
        float edge1 = Mathf.Abs(Vector3.Dot(triangle.Vertex1 - triangle.Vertex2, triangle.Vertex1 - point) / ((triangle.Vertex1 - triangle.Vertex2).magnitude * (triangle.Vertex1 - point).magnitude));
        float edge2 = Mathf.Abs(Vector3.Dot(triangle.Vertex2 - triangle.Vertex3, triangle.Vertex2 - point) / ((triangle.Vertex2 - triangle.Vertex3).magnitude * (triangle.Vertex2 - point).magnitude));
        float edge3 = Mathf.Abs(Vector3.Dot(triangle.Vertex3 - triangle.Vertex1, triangle.Vertex3 - point) / ((triangle.Vertex3 - triangle.Vertex1).magnitude * (triangle.Vertex3 - point).magnitude));
        float[] edges = new float[] { edge1, edge2, edge3 };

        int edgeIndex = edges.ToList().IndexOf(edges.Max());

        Vector3 waterNormal = Vector3.Cross(waterTriangle.Vertex1 - waterTriangle.Vertex2, waterTriangle.Vertex1 - waterTriangle.Vertex3);

        if (edgeIndex == 0)
        {
            return Vector3.Dot(waterNormal, triangle.Vertex2 - triangle.Vertex1) > 0 ? triangle.Vertex1 : triangle.Vertex2;
        }
        if (edgeIndex == 1)
        {
            return Vector3.Dot(waterNormal, triangle.Vertex3 - triangle.Vertex2) > 0 ? triangle.Vertex2 : triangle.Vertex3;
        }
        else
        {
            return Vector3.Dot(waterNormal, triangle.Vertex1 - triangle.Vertex3) > 0 ? triangle.Vertex3 : triangle.Vertex1;
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

    private float TriangleHeight(Triangle triangle) //Find vertical position of triangle
    {
        return PointHeight(transform.TransformPoint((triangle.Vertex1 + triangle.Vertex2 + triangle.Vertex3) / 3f));
    }
    private float MaxTriangleHeight(Triangle triangle)
    {
        return Mathf.Max(Mathf.Max(
                        PointHeight(transform.TransformPoint(triangle.Vertex1)),
                        PointHeight(transform.TransformPoint(triangle.Vertex2))),
                        PointHeight(transform.TransformPoint(triangle.Vertex3)));
    }
    private (float, Vector3) MinTriangleHeight(Triangle triangle)
    {
        List<float> heightValues = new List<float>();
        heightValues.Add(PointHeight(transform.TransformPoint(triangle.Vertex1)));
        heightValues.Add(PointHeight(transform.TransformPoint(triangle.Vertex2)));
        heightValues.Add(PointHeight(transform.TransformPoint(triangle.Vertex3)));

        float minHeight = heightValues.Min();
        int index = heightValues.IndexOf(minHeight);

        return (minHeight, triangle.GetVertices()[index]);
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


    private Triangle[] TriangulateIntermediate()
    {
        throw new NotImplementedException();
    }

    #endregion

    #region Compute forces on rigidbody ==============================================================================

    private void ApplyBuoyancy(Triangle[] bottomHalf, Triangle[] trianglesSea)
    {
        (Vector3 barycentreBoat, float volumeBoat, Vector3 linearDragBoat, Vector3 angularDragBoat) = MeshHelper.ComputeVolumeAndBarycentre(bottomHalf, transform.position, transform.localToWorldMatrix, rigidbody.velocity, rigidbody.angularVelocity, C);
        (Vector3 barycentreSea, float volumeSea, Vector3 linearDragSea, Vector3 angularDragSea) = MeshHelper.ComputeVolumeAndBarycentre(trianglesSea, transform.position, Matrix4x4.identity, rigidbody.velocity, rigidbody.angularVelocity, C, 1.2f, 1E-5f);

        newBarycentre = (barycentreBoat * volumeBoat + barycentreSea * volumeSea) / (volumeBoat + volumeSea);
        newVolume = volumeSea + volumeBoat;

        rigidbody.AddForceAtPosition((volumeSea + volumeBoat) * 1000 * 9.81f * waterIntersectionNormal, newBarycentre);
        rigidbody.AddForce(linearDragBoat + linearDragSea);

        rigidbody.AddTorque(angularDragSea + angularDragBoat);


        // UnityEngine.Debug.Log("Angular drag boat " + (angularDragBoat).magnitude);
        // UnityEngine.Debug.Log("Angular drag sea " + (angularDragSea).magnitude);
        // UnityEngine.Debug.Log("Angular speed " + rigidbody.angularVelocity.magnitude);
        // UnityEngine.Debug.DrawRay(transform.position, angularDragBoat, Color.blue);
        // UnityEngine.Debug.DrawRay(transform.position, rigidbody.angularVelocity, Color.red);
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
                Gizmos.DrawSphere(chain[i + 1], 0.01f);
                Gizmos.DrawLine(chain[i], chain[i + 1]);
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
            Gizmos.DrawSphere(newBarycentre, 0.1f);
        }

        if (showBelowWaterPoints && meshFilter)
        {
            Gizmos.color = Color.red;
            Gizmos.DrawSphere(averageWaterPosition, 0.1f);
            foreach (Vector3 vertex in meshFilter.mesh.vertices)
            {
                if (PointHeight(transform.TransformPoint(vertex)) < PointHeight(averageWaterPosition))
                {
                    Gizmos.color = Color.blue;
                    Gizmos.DrawSphere(transform.TransformPoint(vertex), 0.01f);
                }
            }
        }

        if (showCenterOfMass && rigidbody)
        {
            Gizmos.color = Color.green;
            Gizmos.DrawSphere(transform.TransformPoint(rigidbody.centerOfMass), 0.1f);
        }

    }

    void DrawTriangleHelper(int i, Mesh mesh, Color color)
    {
        Vector3 v1 = mesh.vertices[mesh.triangles[i * 3]];
        Vector3 v2 = mesh.vertices[mesh.triangles[i * 3 + 1]];
        Vector3 v3 = mesh.vertices[mesh.triangles[i * 3 + 2]];

        Vector3 start = mesh == meshFilter.mesh ? Vector3.zero : Vector3.down * 2;
        UnityEngine.Debug.DrawLine(start, (v1 + v2 + v3) / 3f, color, 5000f);
    }

    void DrawTriangleHelper(Triangle triangle, Color color)
    {

        Vector3 start = Vector3.zero;
        UnityEngine.Debug.DrawLine(start, (triangle.Vertex1 + triangle.Vertex2 + triangle.Vertex3) / 3f, color, 5000f);
    }


    #endregion


}
