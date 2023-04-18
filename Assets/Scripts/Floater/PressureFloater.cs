using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

[RequireComponent(typeof(Rigidbody), typeof(MeshCollider))]
public class PressureFloater : MonoBehaviour
{
    #region Variables =====================================================================================

    bool meshReady; //Set after the event

    //Data structures and vectors (global for optimization purposes)
    Triangle[] waterTriangleNeighbors;
    Cell[] gridCells;
    Triangle[] boatTriangleNeighbors;
    Vertex[] boatVertices;

    List<Triangle> triangleCandidateWater;


    int kernelTriangleCandidate;

    ComputeBuffer tempBuffer;

    MovingAverage smoothedVector;
    Vector3[] waterVertices;
    int[] waterTriangles;

    private float boatMaxSize;

    //References ================
    Rigidbody boatRigidbody;
    WaterManager water;


    //Editor variables ================================
    [Header("References")]
    ComputeShader triangleCandidateShader;

    [SerializeField]
    private Mesh floatingMesh;

    [Header("Parameters")]
    public Vector3 cellMaxSize = new Vector3(0.2f, 0.2f, 0.2f);
    public Vector3 centerOfMass; // To adjust the center of mass if wanted
    public float C = 1;
    public int waterGridWidth = 5;
    public bool useColliderAsMesh = true; // Wether to get the mesh from the collider and use it as the floating mesh

    [Header("Debug option")]
    public bool showCandidatesBoat;
    public bool showCandidatesWater;
    public bool showDragForces;
    public bool showUnderwaterMesh;
    public bool showBouyancyForce;
    public bool showVerticesUnderwater;

    #endregion


    void Start()
    {
        water = GameObject.FindObjectOfType<WaterManager>();
        triangleCandidateShader = StaticResourcesLoader.TriangleCandidateShader;
        triangleCandidateShader = Instantiate(triangleCandidateShader);

        boatRigidbody = GetComponent<Rigidbody>();
        boatRigidbody.centerOfMass = centerOfMass;

        MeshCollider meshCollider = GetComponent<MeshCollider>();
        meshCollider.convex = true;

        if (useColliderAsMesh)
        {
            floatingMesh = meshCollider.sharedMesh;
        }

        triangleCandidateWater = new List<Triangle>();

        if (floatingMesh != null) MeshDataPrecomputation();
    }


    void FixedUpdate()
    {
        doneTrianglePairs.Clear();
        doneTriangleCandidate.Clear();

        if (!meshReady) return;
        UpdateWaterMesh();


        HashSet<Triangle> triangleCandidateBoat = FindTriangleCandidatesGPU();

        // DebugHelper.ShowMesh(boatTriangleNeighbors, transform, Color.blue);
        // DebugHelper.ShowMesh(triangleCandidateBoat.ToArray(), transform, Color.blue);

        if (showCandidatesBoat) DebugHelper.ShowMesh(triangleCandidateBoat.ToArray(), transform, Color.red);
        if (showCandidatesWater) DebugHelper.ShowMesh(triangleCandidateWater.ToArray(), transform, Color.cyan);

        ComputeVerticesDepth();

        ApplyBuoyancy(triangleCandidateBoat, boatTriangleNeighbors);
    }

    #region Precomputations ======================================================================================
    private void MeshDataPrecomputation() //Executing once the mesh has been created
    {
        floatingMesh = MeshHelper.WeldVertices(floatingMesh);
        // (Vector3 barycentre, float fullVolume) = MeshHelper.ComputeVolumeAndBarycentre(floatingMesh.vertices, floatingMesh.triangles, transform);

        //Precomputing the triangle neighbor relations and the background grid (only once at the beginning because we suppose that the mesh does not change)
        (boatTriangleNeighbors, boatVertices) = MeshHelper.FindTriangleNeighbors(floatingMesh);

        gridCells = MeshHelper.ComputeBackgroundGrid(floatingMesh.bounds, cellMaxSize, floatingMesh.triangles.Length);

        AssignTrianglesToBoundingBox();

        CreateWaterGrid();
        (waterTriangleNeighbors, _) = MeshHelper.FindTriangleNeighbors(waterVertices, waterTriangles);

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
        Vector3 mins, maxes;
        Bounds triangleBoundingbox;

        foreach (Triangle triangle in boatTriangleNeighbors)
        {
            Vector3[] positions = triangle.GetVerticesPosition();

            mins = Vector3.Min(Vector3.Min(positions[0], positions[1]), positions[2]);
            maxes = Vector3.Max(Vector3.Max(positions[0], positions[1]), positions[2]);

            triangleBoundingbox = new Bounds((maxes + mins) / 2, maxes - mins);

            for (int j = 0; j < gridCells.Length; j++)
            {
                if (gridCells[j].Intersects(triangleBoundingbox))
                {
                    gridCells[j].AddSet1(triangle);
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



    private void UpdateWaterMesh()
    {
        for (int i = 0, y = 0; y <= waterGridWidth; y++)
        {
            for (int x = 0; x <= waterGridWidth; x++, i++)
            {
                waterVertices[i] = transform.TransformPoint(new Vector3(x / (float)waterGridWidth * boatMaxSize - boatMaxSize / 2, 0, y / (float)waterGridWidth * boatMaxSize - boatMaxSize / 2)); //World coordinates
                waterVertices[i] = transform.InverseTransformPoint(new Vector3(waterVertices[i].x, water.GetWaterHeight(waterVertices[i]), waterVertices[i].z));
            }
        }

        for (int i = 0; i < waterTriangleNeighbors.Length; i++)
        {
            waterTriangleNeighbors[i].Vertex1Pos = waterVertices[waterTriangles[waterTriangleNeighbors[i].N1]];
            waterTriangleNeighbors[i].Vertex2Pos = waterVertices[waterTriangles[waterTriangleNeighbors[i].N2]];
            waterTriangleNeighbors[i].Vertex3Pos = waterVertices[waterTriangles[waterTriangleNeighbors[i].N3]];
        }
    }

    private void ComputeVerticesDepth()
    {
        for (int i = 0; i < boatVertices.Length; i++)
        {
            boatVertices[i].depth = transform.TransformPoint(boatVertices[i].position).y - water.GetWaterHeight(boatVertices[i].position);
        }
    }

    private HashSet<Triangle> FindTriangleCandidatesGPU() //Optimized version of the cpu one above
    {
        //Set the frame data
        int triangleCount = waterTriangles.Length / 3;
        tempBuffer = ComputeHelper.CreateAndSetBuffer(waterVertices, triangleCandidateShader, "vertices", kernelTriangleCandidate);
        ComputeBuffer bufferCombined = ComputeHelper.CreateAndSetBuffer<int>((triangleCount + 1) * gridCells.Length, triangleCandidateShader, "combinedBuffer", kernelTriangleCandidate);

        ComputeHelper.Run(triangleCandidateShader, gridCells.Length, 1, 1, kernelTriangleCandidate);

        int[] combined = new int[(triangleCount + 1) * gridCells.Length];
        bufferCombined.GetData(combined); // Freezing line
        bufferCombined.Release();

        tempBuffer.Release();

        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>();
        if (showCandidatesWater) triangleCandidateWater.Clear();

        for (int i = 0; i < gridCells.Length; i++)
        {
            gridCells[i].resetSet2();
            if (combined[i * (triangleCount + 1)] > 0)
            {
                triangleCandidateBoatHS.UnionWith(gridCells[i].TriangleSet1);
                for (int j = 0; j < combined[i * (triangleCount + 1)]; j++)
                {
                    Triangle triangle = waterTriangleNeighbors[combined[i * (triangleCount + 1) + j + 1]];
                    gridCells[i].AddSet2(triangle);
                    if (showCandidatesWater) triangleCandidateWater.Add(triangle);
                }
            }
        }

        return triangleCandidateBoatHS;
    }

    #region Compute forces on rigidbody ==============================================================================
    private HashSet<(Triangle, Triangle)> doneTrianglePairs = new HashSet<(Triangle, Triangle)>();
    private HashSet<Triangle> doneTriangleCandidate = new HashSet<Triangle>();

    private void ApplyBuoyancy(HashSet<Triangle> triangleCandidate, Triangle[] trianglesBoat)
    {
        List<Triangle> underwaterTriangles = new List<Triangle>();

        foreach (Triangle triangle in trianglesBoat)
        {
            if (triangle.GetVertices()[0].depth > 0 || triangle.GetVertices()[1].depth > 0 || triangle.GetVertices()[2].depth > 0) continue;
            if (triangleCandidate.Contains(triangle)) continue;

            underwaterTriangles.Add(triangle);

        }

        foreach (Cell cell in gridCells)
        {
            if (!cell.HasCandidates()) continue;
            // DebugHelper.ShowMesh(cell.TriangleSet1, transform, Color.cyan);
            foreach (Triangle triangleSet2 in cell.TriangleSet2)
            {
                foreach (Triangle triangleSet1 in cell.TriangleSet1)
                {
                    if (doneTrianglePairs.Contains((triangleSet1, triangleSet2))) continue;
                    doneTrianglePairs.Add((triangleSet1, triangleSet2));

                    Triangle[] splitTriangles = GetSplitTriangle(triangleSet1, triangleSet2);


                    foreach (Triangle triangle in splitTriangles)
                    {
                        underwaterTriangles.Add(triangle);
                    }
                }
            }
        }


        Vector3 angularVelocity = transform.InverseTransformDirection(boatRigidbody.angularVelocity);
        Vector3 linearVelocity = transform.InverseTransformDirection(boatRigidbody.velocity);
        (Vector3 barycentreBoat, float volumeBoat, List<Vector3> forcePositions, List<Vector3> forceDirections) = MeshHelper.ComputeVolumeAndBarycentre(underwaterTriangles.ToArray(), linearVelocity, angularVelocity, C);

        boatRigidbody.AddForceAtPosition(volumeBoat * 1000 * 9.81f * Vector3.up, transform.TransformPoint(barycentreBoat));

        for (int i = 0; i < forcePositions.Count; i++)
        {
            boatRigidbody.AddForceAtPosition(transform.TransformDirection(forceDirections[i]), transform.TransformPoint(forcePositions[i]));
        }

        if (showUnderwaterMesh) DebugHelper.ShowMesh(underwaterTriangles.ToArray(), transform, Color.green);
        if (showBouyancyForce) Debug.DrawRay(transform.TransformPoint(barycentreBoat), volumeBoat * 1000 * 9.81f * Vector3.up, Color.red);
        if (showDragForces)
        {
            Vector3 overallForce = Vector3.zero;
            float maxLength = Mathf.Max(forceDirections.Select(i => i.magnitude).ToArray());

            for (int i = 0; i < forcePositions.Count; i++)
            {
                Color color = Color.Lerp(Color.black, Color.red, forceDirections[i].magnitude / maxLength);
                Debug.DrawRay(transform.TransformPoint(forcePositions[i]), transform.TransformDirection(forceDirections[i] * 0.01f), color); // + new Vector3(Random.Range(-0.01f, 0.01f), Random.Range(-0.01f, 0.01f), Random.Range(-0.01f, 0.01f))
                overallForce += transform.TransformDirection(forceDirections[i]);
            }

            Debug.DrawRay(transform.position, overallForce, Color.cyan);
        }
    }


    private (List<Vector3>, List<Vector3>) ComputeForceAtTriangle(Triangle triangle)
    {
        List<Vector3> forcePositions = new List<Vector3>();
        List<Vector3> forceDirections = new List<Vector3>();

        float rho = 1000f;
        float mu = 1E-3f;
        int s = 4;

        Vertex[] vertices = triangle.GetVertices();
        Vector3[] verticesPosition = triangle.GetVerticesPosition();
        Vector3 normal = Vector3.Cross(verticesPosition[0] - verticesPosition[1], verticesPosition[0] - verticesPosition[2]);

        Vector3 angularSpeed = transform.InverseTransformDirection(boatRigidbody.angularVelocity);
        Vector3 linearSpeed = transform.InverseTransformDirection(boatRigidbody.velocity);

        for (int i = 0; i < s; i++)
        {
            for (int j = 0; j < s - i; j++)
            {
                float a = (i + 0.25f) / (float)s;
                float b = (j + 0.25f) / (float)s;
                Vector3 position = verticesPosition[0] + a * (verticesPosition[1] - verticesPosition[0]) + b * (verticesPosition[2] - verticesPosition[0]);
                float depth = vertices[0].depth + a * (vertices[1].depth - vertices[0].depth) + b * (vertices[2].depth - vertices[0].depth);

                float pressure = rho * 9.81f * depth;
                float surfaceArea = normal.magnitude * 0.5f / (s * (s + 1) * 0.5f);

                Vector3 forceOnTriangle = normal.normalized * pressure * surfaceArea;
                Vector3 speedAtTriangle = Vector3.Cross(angularSpeed, position);


                //Friction torque
                float projectedSurfaceAngularSpeed = Vector3.Dot(normal, angularSpeed) > 0 ? ProjectedSurface(verticesPosition[0], verticesPosition[1], verticesPosition[2], speedAtTriangle) / (s * (s + 1) * 0.5f) : 0;
                Vector3 angularDrag = -rho * 0.5f * speedAtTriangle.magnitude * speedAtTriangle.magnitude *
                                projectedSurfaceAngularSpeed * position.magnitude * angularSpeed.normalized * C;

                //Shear stress
                float projectedSurfaceShearStress = ProjectedSurface(verticesPosition[0], verticesPosition[1], verticesPosition[2], position) / (s * (s + 1) * 0.5f);
                angularDrag -= projectedSurfaceShearStress * mu * angularSpeed * position.magnitude;

                //Linear friction drag
                float projectedSurfaceLinearSpeed = Vector3.Dot(normal, linearSpeed) > 0 ? ProjectedSurface(verticesPosition[0], verticesPosition[1], verticesPosition[2], linearSpeed) / (s * (s + 1) * 0.5f) : 0;
                Vector3 linearDrag = -rho * 0.5f * linearSpeed.magnitude * linearSpeed * projectedSurfaceLinearSpeed * C;

                forcePositions.Add(position);
                forceDirections.Add(forceOnTriangle);
            }
        }

        return (forcePositions, forceDirections);
    }

    private static float ProjectedSurface(Vector3 v1, Vector3 v2, Vector3 v3, Vector3 normal)
    {
        Vector3 vertex1Projected = v1 - Vector3.Dot(v1, normal.normalized) * normal.normalized;
        Vector3 vertex2Projected = v2 - Vector3.Dot(v2, normal.normalized) * normal.normalized;
        Vector3 vertex3Projected = v3 - Vector3.Dot(v3, normal.normalized) * normal.normalized;
        return (Vector3.Cross(vertex1Projected - vertex2Projected, vertex1Projected - vertex3Projected)).magnitude / 2f;
    }

    #endregion

    #region Split triangles ===================================================================

    private Triangle[] GetSplitTriangle(Triangle triangleSet1, Triangle triangleSet2)
    {
        List<Triangle> splitTriangles = new List<Triangle>();


        Vector3 a = triangleSet1.Vertex1Pos; // a, b, c is the boat
        Vector3 b = triangleSet1.Vertex2Pos;
        Vector3 c = triangleSet1.Vertex3Pos;
        Vector3 d = triangleSet2.Vertex1Pos;
        Vector3 e = triangleSet2.Vertex2Pos;
        Vector3 f = triangleSet2.Vertex3Pos;


        List<Vector3?> intersectionsWater = new List<Vector3?>();
        List<Vector3?> intersectionsBoat = new List<Vector3?>();

        intersectionsWater.Add(IntersectionEdgeTriangle(a, b, c, d, e));
        intersectionsWater.Add(IntersectionEdgeTriangle(a, b, c, e, f));
        intersectionsWater.Add(IntersectionEdgeTriangle(a, b, c, f, d));
        int intersectionPointsCounterWater = intersectionsWater.Count(s => s != null);

        intersectionsBoat.Add(IntersectionEdgeTriangle(d, e, f, a, b));
        intersectionsBoat.Add(IntersectionEdgeTriangle(d, e, f, b, c));
        intersectionsBoat.Add(IntersectionEdgeTriangle(d, e, f, c, a));
        int intersectionPointsCounterBoat = intersectionsBoat.Count(s => s != null);

        // Debug.Log(intersectionPoints.Count);
        Vertex[] vertices = triangleSet1.GetVertices();

        if (intersectionPointsCounterBoat + intersectionPointsCounterWater == 0 && vertices[0].depth < 0 && vertices[1].depth < 0 && vertices[2].depth < 0 && !doneTriangleCandidate.Contains(triangleSet1))
        {
            splitTriangles.Add(triangleSet1);
            doneTriangleCandidate.Add(triangleSet1);
        }
        else if (intersectionPointsCounterBoat + intersectionPointsCounterWater == 2)
        {
            List<Vector3> tempList = new List<Vector3>();
            for (int i = 0; i < 3; i++)
            {
                if (intersectionsBoat[i] is not null) tempList.Add((Vector3)intersectionsBoat[i]);
                if (intersectionsWater[i] is not null) tempList.Add((Vector3)intersectionsWater[i]);
            }
            // Debug.DrawLine(transform.TransformPoint(tempList[0]), transform.TransformPoint(tempList[1]), Color.red);

            // If we have two intersection points, it means we have a proper split
            List<Vertex> polygon = new List<Vertex>();
            // Case 1: the water triangle in larger than the boat triangle (the two intersection points are on edges of the boat triangle)
            if (intersectionPointsCounterBoat == 2)
            {
                for (int i = 0; i < vertices.Length; i++)
                {
                    if (vertices[i].depth <= 0) polygon.Add(vertices[i]);
                    if (intersectionsBoat[i] is not null) polygon.Add(new Vertex((Vector3)intersectionsBoat[i], 0));
                }
            }
            // Case 2: the boat triangle in larger than the water triangle (the two intersection points are inside of the boat triangle)
            else if (intersectionPointsCounterWater == 2)
            {
                List<Vertex> underwaterVertices = vertices.Where(item => item.depth <= 0).ToList();
                if (underwaterVertices.Count > 0)
                {
                    foreach (Vector3? vertex in intersectionsWater.Where(s => s is not null))
                    {
                        polygon.Add(new Vertex((Vector3)vertex, 0));
                    }
                    Vertex averageVertex = underwaterVertices.Aggregate((v1, v2) => new Vertex(v1.position + v2.position, 0));
                    averageVertex.position = averageVertex.position / underwaterVertices.Count();
                    polygon.Add(averageVertex);
                }

            }
            // Case 3: one of the intersection point is on an edge of the boat triangle and the other one is inside the boat triangle
            else
            {
                polygon.Add(new Vertex((Vector3)intersectionsBoat.Where(s => s is not null).ToList()[0], 0));

                for (int i = 0; i < 3; i++)
                {
                    if (intersectionsBoat[i] is null) continue;

                    if (vertices[i].depth <= 0) polygon.Add(vertices[i]);
                    else polygon.Add(vertices[(i + 1) % vertices.Length]);
                    break;
                }

                List<Vertex> underwaterVertices = vertices.Where(item => item.depth <= 0).ToList();
                if (underwaterVertices.Count > 0)
                {
                    Vertex averageVertex = underwaterVertices.Aggregate((v1, v2) => new Vertex(v1.position + v2.position, 0));
                    averageVertex.position = averageVertex.position / underwaterVertices.Count();
                    if (averageVertex.position != polygon[polygon.Count - 1].position) polygon.Add(averageVertex);

                    polygon.Add(new Vertex((Vector3)intersectionsWater.Where(s => s is not null).ToList()[0], 0));
                }
            }

            // Triangulate the polygon (one or two triangles depending on the shape)
            for (int i = 0; i < polygon.Count; i += 2)
            {
                if (Vector3.Dot(Vector3.Cross(polygon[i].position - polygon[(i + 1) % polygon.Count].position, polygon[i].position - polygon[(i + 2) % polygon.Count].position), triangleSet1.GetNormal()) > 0)
                {
                    splitTriangles.Add(new Triangle(0, polygon[i], polygon[(i + 1) % polygon.Count], polygon[(i + 2) % polygon.Count]));
                }
                else splitTriangles.Add(new Triangle(0, polygon[(i + 1) % polygon.Count], polygon[i], polygon[(i + 2) % polygon.Count]));
            }

        }
        // Should never happen (it seems it doesn't which is good)
        else if (intersectionPointsCounterBoat + intersectionPointsCounterWater == 3 || intersectionPointsCounterBoat + intersectionPointsCounterWater == 1)
        {
            Debug.Log(intersectionPointsCounterBoat);
            Debug.Log(intersectionPointsCounterWater + "\n");
        }

        return splitTriangles.ToArray();
    }

    public static Vector3? IntersectionEdgeTriangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 e)
    {
        Vector3 n = Vector3.Cross(a - b, a - c).normalized;

        float E = Vector3.Dot(a - e, n);
        float D = Vector3.Dot(a - d, n);
        //Are d and e on different side of triangle abc ?
        if (D * E > 0) return null;

        float t = D / (D - E);
        Vector3 p = t * e + (1 - t) * d;

        // if ((a - b).normalized == (a - p).normalized ||
        //     (b - c).normalized == (b - p).normalized ||
        //     (c - a).normalized == (c - p).normalized)
        // {
        //     Debug.Log("On edge");
        // }

        if (((a - b).normalized != (a - p).normalized && Vector3.Dot(Vector3.Cross((a - b).normalized, (a - p).normalized), n) < 0) ||
            ((b - c).normalized != (b - p).normalized && Vector3.Dot(Vector3.Cross((b - c).normalized, (b - p).normalized), n) < 0) ||
            ((c - a).normalized != (c - p).normalized && Vector3.Dot(Vector3.Cross((c - a).normalized, (c - p).normalized), n) < 0)) return null;

        return p;
    }

    #endregion

    public void OnDrawGizmos()
    {
        if (!Application.isPlaying || !showVerticesUnderwater) return;

        foreach (Vertex vertex in boatVertices)
        {
            Gizmos.color = vertex.depth > 0 ? Color.red : Color.blue;
            Gizmos.DrawSphere(transform.TransformPoint(vertex.position), 0.05f);


        }
    }
}

