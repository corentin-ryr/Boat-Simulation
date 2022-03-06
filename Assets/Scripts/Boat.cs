using UnityEngine;
using UnityEditor;
using System.Collections.Generic;
using System;
using System.Collections;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using System.Diagnostics;
using System.Linq;
using Habrador_Computational_Geometry;

[RequireComponent(typeof(Rigidbody))]
public class Boat : MonoBehaviour
{
    #region Variables =====================================================================================

    //Internal variables of the mesh
    float volume;
    Vector3 barycentre;
    Vector3 newBarycentre;

    bool meshReady; //Set after the event

    //Data structures and vectors (global for optimization purposes)
    Triangle[] boatTriangleNeighbors;
    Triangle[] waterTriangleNeighbors;
    Cell[] gridCells;
    Triangle[] triangleCandidateBoat;
    Triangle[] triangleCandidateWater;
    Vector3[] projectedVertices;

    List<Vector3> chain;
    List<Vector3> lowerChain;
    Vector3 a, b, c, d, e, f;
    Vector3 waterIntersectionNormal;
    private Vector3 averageWaterPosition;

    //References
    MeshFilter meshFilter;
    Rigidbody rigidbody;
    Stopwatch stopWatch;

    //Editor variables
    [Header("References")]
    public Water water;

    [Header("Parameters")]
    public Vector3 cellMaxSize = new Vector3(0.2f, 0.2f, 0.2f);
    public bool procedurallyGeneratedMesh;
    public Vector3 centerOfMass;



    [Header("Debug option")]
    public bool showBackgroundGrid;
    public bool showStatValues;
    public bool showCandidatesBoat;
    public bool showCandidatesWater;
    public bool showIntersectionLine;
    public bool showIntersectionMesh;
    public bool showLowerChain;
    public bool showWaterNormal;



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
        stopWatch = new Stopwatch();

        chain = new List<Vector3>();
        lowerChain = new List<Vector3>();

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
            UnityEngine.Debug.Log(volume);
            UnityEngine.Debug.Log(barycentre);
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

        stopWatch.Restart();
        int[] cellCandidates = FindTriangleCandidates();


        if (showCandidatesBoat)
        {
            DebugHelper.ShowMesh(triangleCandidateBoat, transform, Color.red);
        }
        if (showCandidatesWater)
        {
            DebugHelper.ShowMesh(triangleCandidateWater, water.transform, Color.blue);
        }

        stopWatch.Stop();
        UnityEngine.Debug.Log("Time for finding candidates " + stopWatch.Elapsed.Milliseconds);


        stopWatch.Restart();

        //TODO sometimes it bugs
        ComputeIntersectingLine(cellCandidates);


        if (chain?.Count == 0)
        {
            //Test if there is water above
            UnityEngine.Debug.DrawRay(transform.position, Vector3.up * 10);
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

        stopWatch.Stop();
        UnityEngine.Debug.Log("Time to find intersection line & bottom half " + stopWatch.Elapsed.Milliseconds);
        stopWatch.Restart();


        if (bottomHalf.Length == 0) return;

        //We triangulate the sea (we just create a plane with the circle)
        List<Triangle> seaTriangulated = new List<Triangle>();
        for (int i = 0; i < chain.Count-1; i++)
        {
            seaTriangulated.Add(new Triangle(0, averageWaterPosition, chain[i], chain[i + 1]));
        }


        if (showIntersectionMesh)
        {
            DebugHelper.ShowMesh(bottomHalf, transform, Color.blue);
            DebugHelper.ShowMesh(seaTriangulated.ToArray(), Color.blue);
        }


        stopWatch.Stop();
        UnityEngine.Debug.Log("Time for recontructing mesh " + stopWatch.Elapsed.Milliseconds);

        //Apply the forces
        ApplyBuoyancy(bottomHalf, seaTriangulated.ToArray());
    }


    #region Precomputations ======================================================================================
    private void MeshDataPrecomputation() //Executing once the mesh has been created
    {
        meshFilter.mesh = MeshHelper.WeldVertices(meshFilter.mesh, 1E-5f);
        (barycentre, volume) = MeshHelper.ComputeVolumeAndBarycentre(meshFilter.mesh.vertices, meshFilter.mesh.triangles);

        //Precomputing the triangle neighbor relations and the background grid (only once at the beginning because we suppose that the mesh does not change)
        boatTriangleNeighbors = MeshHelper.FindTriangleNeighbors(meshFilter.mesh);

        gridCells = MeshHelper.ComputeBackgroundGrid(meshFilter.mesh.bounds, cellMaxSize, meshFilter.mesh.triangles.Length);

        AssignTrianglesToBoundingBox();

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
    }
    #endregion

    #region Intersection computation with water ======================================================

    private int[] FindTriangleCandidates() //Return the indices of the cells with candidate triangle (for finding the first intersection of the chain)
    {
        //Assign water triangles to boat bounding boxes
        Vector3[] vertices = water.Mesh.vertices;
        int[] triangles = water.Mesh.triangles;
        HashSet<int> candidateCellsHS = new HashSet<int>();

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
                    candidateCellsHS.Add(j);
                }
            }
        }

        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>(meshFilter.mesh.triangles.Length / 3);
        HashSet<Triangle> triangleCandidateWaterHS = new HashSet<Triangle>(triangles.Length / 3);

        int[] candidateCells = new int[candidateCellsHS.Count];
        candidateCellsHS.CopyTo(candidateCells);
        foreach (int cellIndex in candidateCells)
        {
            triangleCandidateBoatHS.UnionWith(gridCells[cellIndex].TriangleSet1);
            triangleCandidateWaterHS.UnionWith(gridCells[cellIndex].TriangleSet2);
        }

        triangleCandidateBoat = new Triangle[triangleCandidateBoatHS.Count];
        triangleCandidateBoatHS.CopyTo(triangleCandidateBoat);
        triangleCandidateWater = new Triangle[triangleCandidateWaterHS.Count];
        triangleCandidateWaterHS.CopyTo(triangleCandidateWater);

        return candidateCells;
    }

    private void ComputeIntersectingLine(int[] cellCandidates)
    {
        chain.Clear();
        // lowerChain.Clear();

        //Find first intersection
        bool swaped = false;
        Triangle currentTriangle;
        Triangle triangleToCheck;
        Vector3 currentP = Vector3.positiveInfinity;
        Triangle[] nextCurrentTriangles = new Triangle[0];

        foreach (int cellIndex in cellCandidates)
        {
            Cell cell = gridCells[cellIndex];
            for (int i = 0; i < cell.TriangleSet2.Length; i++) //The set of triangles in the cell (water triangles)
            {
                triangleToCheck = triangleCandidateWater[i];
                for (int j = 0; j < cell.TriangleSet1.Length; j++) //The set of triangles in the cell (boat triangles)
                {
                    currentTriangle = triangleCandidateBoat[j];
                    bool worked = false;
                    (worked, nextCurrentTriangles) = IdentifyCase(currentTriangle, meshFilter.mesh, triangleToCheck, water.Mesh, ref currentP, ref swaped, ref triangleToCheck);
                    if (worked)
                    {
                        chain.Add(currentP);
                        goto End;
                    }
                }
            }
        }

        return;
    End:

        int tempI = 0;
        bool loopContiue = true;
        //We have the first triangle, we start looping
        while (loopContiue)
        {
            foreach (Triangle nextCurrentTriangle in nextCurrentTriangles) //Normally only one triangle in this array but if we have cases 3 or 4 we can have 0 or more than 1 triangles.
            {
                currentTriangle = nextCurrentTriangle;
                (loopContiue, nextCurrentTriangles) = IdentifyCase(currentTriangle, swaped ? water.Mesh : meshFilter.mesh, triangleToCheck,
                                                                    swaped ? meshFilter.mesh : water.Mesh, ref currentP, ref swaped, ref triangleToCheck);
                if (loopContiue)
                {
                    chain.Add(currentP);
                    break;
                }
                UnityEngine.Debug.Log("No intersection");
                return;
            }

            tempI++;
            if ((currentP - chain[0]).magnitude < 1E-3 || tempI > 50) loopContiue = false;
        }
    }

    //Triangle1 is from the boat and triangle2 is from the water
    //TODO compute cases 3 and 4
    private (bool, Triangle[]) IdentifyCase(Triangle triangle1, Mesh mesh1, Triangle triangle2, Mesh mesh2, ref Vector3 P, ref bool swaped, ref Triangle triangleToCheck)
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
            triangleToCheck = swaped ? triangleIndicesToCheck2[0] : triangleIndicesToCheck2[0];
            swaped = !swaped;
            P = tempP1;
            return (true, triangleIndicesToCheck1); //We act like a case 1 but we have a special triangleToCheck
        }
        else if (case1) //Case 1, we have to swap and we set the triangle to check next
        {
            swaped = !swaped;
            triangleToCheck = triangle1;
            P = tempP1;
            return (true, triangleIndicesToCheck1);
        }
        else if (case2)//Case 2, we don't swap and we also set the next triangle to check
        {
            P = tempP2;
            return (true, triangleIndicesToCheck2);
        }


        return (false, triangleIndicesToCheck1);
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


    #endregion

    #region Triangulate the mesh =================================================================================
    private int[] DelaunayTriangulation(Vector3[] verticesToTriangulate)
    {
        projectedVertices = projectPoints(verticesToTriangulate, 1);

        // Vertex is TriangleNet.Geometry.Vertex
        TriangleNet.Geometry.Polygon polygon = new TriangleNet.Geometry.Polygon();
        for (int j = 0; j < projectedVertices.Length; j++)
        {
            polygon.Add(new TriangleNet.Geometry.Vertex(projectedVertices[j].x, projectedVertices[j].z));
        }

        // ConformingDelaunay is false by default; this leads to ugly long polygons at the edges
        // because the algorithm will try to keep the mesh convex
        TriangleNet.Meshing.ConstraintOptions options =
            new TriangleNet.Meshing.ConstraintOptions() { ConformingDelaunay = true, SegmentSplitting = 2 };
        IMesh triangleMesh = (TriangleNet.Mesh)polygon.Triangulate(options);

        List<int> triangleList = new List<int>();
        foreach (TriangleNet.Topology.Triangle triangle in triangleMesh.Triangles)
        {
            triangleList.Add(triangle.vertices[0].id);
            triangleList.Add(triangle.vertices[1].id);
            triangleList.Add(triangle.vertices[2].id);
        }

        //Fill the triangle array of ints used by the unity mesh
        return triangleList.ToArray();

    }

    private Vector3[] projectPoints(Vector3[] vertices, float radius)
    {
        // Choose a point from which to project
        Vector3 projectPoint = transform.up * radius;

        //Elevate the sphere (to be just on the plane)
        Vector3[] projectedVertices = new Vector3[vertices.Length]; // new Vector3[nbPoints];

        //Project the vertices on the plane
        for (int i = 0; i < projectedVertices.Length; i++)
        {
            Vector3 point = vertices[i];
            projectedVertices[i] = new Vector3((point.x * radius) / (radius - point.y), 0, (point.z * radius) / (radius - point.y));
        }

        return projectedVertices;
    }

    private Triangle[] FindBottomHalfOfBoat()
    {
        HashSet<Triangle> bottomHalf = new HashSet<Triangle>();

        bool stillInCandidates = true;
        Triangle currentTriangle = triangleCandidateBoat[0];
        float minDepth = TriangleHeight(currentTriangle);
        while (stillInCandidates)
        {
            float previousDepth = minDepth;
            //Get the neighbor further down
            foreach (Triangle neighbor in currentTriangle.GetNeighbors())
            {
                float depth = TriangleHeight(neighbor);
                if (depth < minDepth)
                {
                    currentTriangle = neighbor;
                    minDepth = depth;
                }
            }

            if (!(new List<Triangle>(triangleCandidateBoat).Contains(currentTriangle)) || previousDepth == minDepth)
            {
                stillInCandidates = false;
            }
        }

        //We get all the triangle by graph traversal (BFS)
        Queue<Triangle> Q = new Queue<Triangle>();
        Q.Enqueue(currentTriangle);
        bottomHalf.Add(currentTriangle);
        float waterSurfaceDepth = PointHeight(averageWaterPosition);

        while (Q.Count > 0)
        {
            currentTriangle = Q.Dequeue();
            foreach (Triangle neighbor in currentTriangle.GetNeighbors())
            {
                if (MaxTriangleHeight(neighbor) < waterSurfaceDepth && bottomHalf.Add(neighbor))
                {
                    Q.Enqueue(neighbor);
                }
            }
        }
        Triangle[] array = new Triangle[bottomHalf.Count];
        bottomHalf.CopyTo(array);
        return array;
    }

    private float TriangleHeight(Triangle triangle)
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
    private float PointHeight(Vector3 point)
    {
        return Vector3.Dot(point, waterIntersectionNormal);
    }

    #endregion

    #region Compute forces on rigidbody ==============================================================================œ

    private void ApplyBuoyancy(Triangle[] bottomHalf, Triangle[] trianglesSea)
    {

        (Vector3 barycentreBoat, float volumeBoat) = MeshHelper.ComputeVolumeAndBarycentre(bottomHalf);
        (Vector3 barycentreSea, float volumeSea) = MeshHelper.ComputeVolumeAndBarycentre(trianglesSea);

        newBarycentre = transform.TransformPoint(barycentreBoat + barycentreSea);

        rigidbody.AddForceAtPosition((volumeSea + volumeBoat) * 1000 * 9.81f * Vector3.up, newBarycentre);
    }


    #endregion


    #region Gizmos and debug =====================================================================================

    void OnDrawGizmos()
    {
        if (showStatValues)
        {
            Gizmos.DrawSphere(transform.TransformPoint(barycentre), 0.1f);
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
        if (lowerChain.Count > 0 && showLowerChain)
        {
            Gizmos.color = Color.black;
            Gizmos.DrawSphere(transform.TransformPoint(lowerChain[0]), 0.01f);
            for (int i = 0; i < lowerChain.Count - 1; i++)
            {
                Gizmos.DrawSphere(transform.TransformPoint(lowerChain[i + 1]), 0.01f);
            }

        }
        // if (projectedVertices != null)
        // {
        //     foreach (Vector3 item in projectedVertices)
        //     {
        //         Gizmos.DrawSphere(transform.TransformPoint(item), 0.01f);
        //     }
        // }

        Gizmos.color = Color.cyan;
        Gizmos.DrawSphere(newBarycentre, 0.1f);
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
