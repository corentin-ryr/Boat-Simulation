
using UnityEngine;
using UnityEditor;
using System.Collections.Generic;
using System;
using System.Collections;
using TriangleNet.Geometry;
using TriangleNet.Meshing;

public class Boat : MonoBehaviour
{
    float volume;
    Vector3 barycentre;
    Triangle[] triangleNeighbors;
    Cell[] gridCells;
    bool meshReady;

    Triangle[] triangleCandidateBoat;
    Triangle[] triangleCandidateWater;

    Triangle[] waterTriangleNeighbors;

    List<Vector3> chain;

    MeshFilter meshFilter;
    Isocahedron isocahedron;

    [Header("References")]
    public ComputeShader computeShader;
    public Water water;

    [Header("Parameters")]
    public Vector3 cellMaxSize = new Vector3(0.2f, 0.2f, 0.2f);
    public int nbSubdivision;


    [Header("Debug option")]
    public bool showBackgroundGrid;
    public bool showStatValues;
    public bool showCandidatesBoat;
    public bool showCandidatesWater;
    public bool showIntersectionLine;

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
        Mesh mesh = new Mesh();

        isocahedron = new Isocahedron(nbSubdivision, 1f, computeShader);
        isocahedron.CreateSphere();
    }

    void Update()
    {
        if (!meshReady) return;

        if (showStatValues) //Used to check if the volume and barycentre computation is working well
        {
            (barycentre, volume) = MeshHelper.ComputeVolumeAndBarycentre(meshFilter.mesh.vertices, meshFilter.mesh.triangles);
            Debug.Log(volume);
            Debug.Log(barycentre);
        }

        if (showBackgroundGrid)
        {
            Debug.Log("Number of cells " + gridCells.Length);
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

        FindTriangleCandidates();

        if (showCandidatesBoat)
        {
            DebugHelper.ShowMesh(triangleCandidateBoat, transform, meshFilter.mesh, Color.red);
        }
        if (showCandidatesWater)
        {
            DebugHelper.ShowMesh(triangleCandidateWater, water.transform, water.Mesh, Color.blue);
        }

        ComputeIntersectingLine();

        //We triangulate the sea (we just create a plane with the circle)
        chain.Add(Vector3.zero);
        Vector3[] verticesToTriangulate = chain.ToArray();
        int[] intersectionSeaTriangles = DelaunayTriangulation(verticesToTriangulate);
        DebugHelper.ShowMesh(verticesToTriangulate, intersectionSeaTriangles, transform, Color.blue);


        //We triangulate the boat (we nned to find the triangles below the line)
        Triangle[] bottomHalf = FindBottomHalfOfBoat();


    }


    #region Precomputations ======================================================================================
    private void MeshDataPrecomputation() //Executing once the mesh has been created
    {
        //Assigning the mesh
        meshFilter.mesh.vertices = isocahedron.vertices.ToArray();
        List<int> triangles = new List<int>();
        foreach (Vector3Int v3 in isocahedron.faces)
        {
            triangles.Add(v3.x);
            triangles.Add(v3.y);
            triangles.Add(v3.z);
        }
        meshFilter.mesh.triangles = triangles.ToArray();

        //Precomputing the triangle neighbor relations and the background grid (only once at the beginning because we suppose that the mesh does not change)
        triangleNeighbors = MeshHelper.FindTriangleNeighbors(meshFilter.mesh);
        gridCells = MeshHelper.ComputeBackgroundGrid(meshFilter.mesh.bounds, cellMaxSize, triangles.Count);

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
                    gridCells[j].AddSet1(triangleNeighbors[i / 3]);
                }
            }
        }
    }
    #endregion

    #region Intersection computation with water ======================================================

    private void FindTriangleCandidates()
    {
        //Assign water triangles to boat bounding boxes
        Mesh waterMesh = water.Mesh;

        foreach (Cell cell in gridCells)
        {
            cell.resetSet2();
        }

        Matrix4x4 transformation = transform.worldToLocalMatrix * water.transform.localToWorldMatrix;
        Vector3 a, b, c;

        Vector3 mins, maxes;

        Bounds triangleBoundingbox;

        for (int i = 0; i < waterMesh.triangles.Length; i += 3)
        {
            a = transformation.MultiplyPoint3x4(waterMesh.vertices[waterMesh.triangles[i]]);
            b = transformation.MultiplyPoint3x4(waterMesh.vertices[waterMesh.triangles[i + 1]]);
            c = transformation.MultiplyPoint3x4(waterMesh.vertices[waterMesh.triangles[i + 2]]);

            mins = Vector3.Min(Vector3.Min(a, b), c);
            maxes = Vector3.Max(Vector3.Max(a, b), c);

            triangleBoundingbox = new Bounds((maxes + mins) / 2, maxes - mins);

            for (int j = 0; j < gridCells.Length; j++)
            {
                if (gridCells[j].HasCandidatePotential && gridCells[j].Intersects(triangleBoundingbox))
                {
                    gridCells[j].AddSet2(waterTriangleNeighbors[i / 3]);
                }
            }
        }

        HashSet<Triangle> triangleCandidateBoatHS = new HashSet<Triangle>();
        HashSet<Triangle> triangleCandidateWaterHS = new HashSet<Triangle>();
        foreach (Cell cell in gridCells)
        {
            if (cell.HasCandidates())
            {
                foreach (Triangle triangle in cell.TriangleSet1)
                {
                    triangleCandidateBoatHS.Add(triangle);
                }

                foreach (Triangle triangle in cell.TriangleSet2)
                {
                    triangleCandidateWaterHS.Add(triangle);
                }
            }
        }

        triangleCandidateBoat = new Triangle[triangleCandidateBoatHS.Count];
        triangleCandidateBoatHS.CopyTo(triangleCandidateBoat);
        triangleCandidateWater = new Triangle[triangleCandidateWaterHS.Count];
        triangleCandidateWaterHS.CopyTo(triangleCandidateWater);
    }

    private void ComputeIntersectingLine()
    {
        chain = new List<Vector3>();

        //Find first intersection
        bool swaped = false;
        Triangle currentTriangle;
        Triangle triangleToCheck;
        Vector3 currentP = Vector3.positiveInfinity;
        int[] nextCurrentTriangles = new int[0];

        for (int i = 0; i < triangleCandidateWater.Length; i++)
        {
            triangleToCheck = triangleCandidateWater[i];
            for (int j = 0; j < triangleCandidateBoat.Length; j++)
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
        return;
    End:


        int tempI = 0;
        bool loopContiue = true;
        //We have the first triangle, we start looping
        while (loopContiue)
        {
            foreach (int triangleIndex in nextCurrentTriangles) //Normally only one triangle in this array but if we have cases 3 or 4 we can have 0 or more than 1 triangles.
            {
                currentTriangle = swaped ? new List<Triangle>(triangleCandidateWater).Find(p => p.TriangleIndex == triangleIndex) : new List<Triangle>(triangleCandidateBoat).Find(p => p.TriangleIndex == triangleIndex);

                // DrawTriangleHelper(currentTriangle.TriangleIndex, swaped ? water.Mesh : meshFilter.mesh, Color.Lerp(Color.green, Color.black, tempI / 6f));

                (loopContiue, nextCurrentTriangles) = IdentifyCase(currentTriangle, swaped ? water.Mesh : meshFilter.mesh, triangleToCheck,
                                                                    swaped ? meshFilter.mesh : water.Mesh, ref currentP, ref swaped, ref triangleToCheck);
                if (loopContiue)
                {
                    chain.Add(currentP);
                    // Debug.Log("Swaped " + tempI + "  " + swaped);
                    // Debug.Log("Step done");

                    break;
                }
                // Debug.DrawLine(Vector3.zero, currentP);
                // DrawTriangleHelper(currentTriangle.TriangleIndex, swaped ? water.Mesh : meshFilter.mesh, Color.green);
                // DrawTriangleHelper(triangleToCheck.TriangleIndex, !swaped ? water.Mesh : meshFilter.mesh, Color.cyan);
                Debug.Log("No intersection");
                return;

            }

            tempI++;

            if ((currentP - chain[0]).magnitude < 1E-3 || tempI > 50) loopContiue = false;
        }


    }

    //Triangle1 is from the boat and triangle2 is from the water
    //TODO compute cases 3 and 4
    private (bool, int[]) IdentifyCase(Triangle triangle1, Mesh mesh1, Triangle triangle2, Mesh mesh2, ref Vector3 P, ref bool swaped, ref Triangle triangleToCheck)
    {
        Matrix4x4 transformation = mesh1 == meshFilter.mesh ? Matrix4x4.identity : transform.worldToLocalMatrix * water.transform.localToWorldMatrix;
        Vector3 a = transformation.MultiplyPoint3x4(mesh1.vertices[mesh1.triangles[triangle1.TriangleIndex * 3]]);
        Vector3 b = transformation.MultiplyPoint3x4(mesh1.vertices[mesh1.triangles[triangle1.TriangleIndex * 3 + 1]]);
        Vector3 c = transformation.MultiplyPoint3x4(mesh1.vertices[mesh1.triangles[triangle1.TriangleIndex * 3 + 2]]);

        transformation = mesh2 == meshFilter.mesh ? Matrix4x4.identity : transform.worldToLocalMatrix * water.transform.localToWorldMatrix;
        Vector3 d = transformation.MultiplyPoint3x4(mesh2.vertices[mesh2.triangles[triangle2.TriangleIndex * 3]]);
        Vector3 e = transformation.MultiplyPoint3x4(mesh2.vertices[mesh2.triangles[triangle2.TriangleIndex * 3 + 1]]);
        Vector3 f = transformation.MultiplyPoint3x4(mesh2.vertices[mesh2.triangles[triangle2.TriangleIndex * 3 + 2]]);


        //Case 1 if p is inside of triangle1 
        int[] triangleIndicesToCheck1 = null;
        Vector3 tempP1 = Vector3.negativeInfinity;
        bool case1 = false;
        if (IntersectionEdgeTriangle(a, b, c, d, e, ref tempP1, P))
        {
            case1 = true;
            triangleIndicesToCheck1 = new int[] { triangle2.N1 };
        }
        else if (IntersectionEdgeTriangle(a, b, c, e, f, ref tempP1, P))
        {
            case1 = true;
            triangleIndicesToCheck1 = new int[] { triangle2.N2 };
        }
        else if (IntersectionEdgeTriangle(a, b, c, f, d, ref tempP1, P))
        {
            case1 = true;
            triangleIndicesToCheck1 = new int[] { triangle2.N3 };
        }


        //Case 2 if p is inside of triangle2
        int[] triangleIndicesToCheck2 = null;
        Vector3 tempP2 = Vector3.negativeInfinity;
        bool case2 = false;
        if (IntersectionEdgeTriangle(d, e, f, a, b, ref tempP2, P))
        {
            case2 = true;
            triangleIndicesToCheck2 = new int[] { triangle1.N1 };
        }
        else if (IntersectionEdgeTriangle(d, e, f, b, c, ref tempP2, P))
        {
            case2 = true;
            triangleIndicesToCheck2 = new int[] { triangle1.N2 };
        }
        else if (IntersectionEdgeTriangle(d, e, f, c, a, ref tempP2, P))
        {
            case2 = true;
            triangleIndicesToCheck2 = new int[] { triangle1.N3 };
        }


        // Debug.Log("P : " + P);
        // Debug.Log(tempP1);
        // Debug.Log(tempP2);
        // Debug.DrawLine(Vector3.zero, tempP1, Color.cyan);
        // Debug.DrawLine(Vector3.zero, tempP2, Color.cyan);

        if (case1 && case2 && tempP1 == tempP2) //Case when we have the point on the edge of the two triangles. The next triangle to check 
        {
            // Debug.Log("Case special");
            triangleToCheck = swaped ? new List<Triangle>(triangleCandidateWater).Find(p => p.TriangleIndex == triangleIndicesToCheck2[0]) : new List<Triangle>(triangleCandidateBoat).Find(p => p.TriangleIndex == triangleIndicesToCheck2[0]);
            swaped = !swaped;
            P = tempP1;
            return (true, triangleIndicesToCheck1); //We act like a case 1 but we have a special triangleToCheck
        }
        else if (case1) //Case 1, we have to swap and we set the triangle to check next
        {
            // Debug.Log("Case 1");

            swaped = !swaped;
            triangleToCheck = triangle1;
            P = tempP1;
            return (true, triangleIndicesToCheck1);
        }
        else if (case2)//Case 2, we don't swap and we also set the next triangle to check
        {
            // Debug.Log("Case 2");

            // Debug.Log(triangleIndicesToCheck2[0]);
            // DrawTriangleHelper(triangleIndicesToCheck2[0], meshFilter.mesh, Color.green);
            P = tempP2;
            return (true, triangleIndicesToCheck2);
        }


        return (false, triangleIndicesToCheck1);
    }

    private bool IntersectionEdgeTriangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 e, ref Vector3 newP, Vector3 previousP)
    {
        Vector3 n = Vector3.Cross(a - b, a - c);

        //Are d and e on different side of triangle abc ?
        if ((Vector3.Dot(a - d, n) * Vector3.Dot(a - e, n) > 0)) return false;


        float t = Vector3.Dot(a - d, n) / (Vector3.Dot(a - d, n) - Vector3.Dot(a - e, n));
        Vector3 p = t * e + (1 - t) * d;

        if (Vector3.Dot(Vector3.Cross(a - b, a - p), n) < -1E-7 ||
            Vector3.Dot(Vector3.Cross(b - c, b - p), n) < -1E-7 ||
            Vector3.Dot(Vector3.Cross(c - a, c - p), n) < -1E-7) return false;

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
        // Vertex is TriangleNet.Geometry.Vertex
        Polygon polygon = new Polygon();
        for (int j = 0; j < verticesToTriangulate.Length; j++)
        {
            polygon.Add(new Vertex(verticesToTriangulate[j].x, verticesToTriangulate[j].z));
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

    private Triangle[] FindBottomHalfOfBoat()
    {
        HashSet<Triangle> bottomHalf = new HashSet<Triangle>();

        bool stillInCandidates = true;
        Triangle currentTriangle = triangleCandidateBoat[0];
        while (stillInCandidates)
        {
            float minDepth = float.MaxValue;
            //Get the neighbor further down
            currentTriangle = TriangleDepth(currentTriangle.N1, meshFilter.mesh) < minDepth ? new List<Triangle>(triangleCandidateBoat).Find(p => p.TriangleIndex == currentTriangle.N1) : currentTriangle;
            currentTriangle = TriangleDepth(currentTriangle.N2, meshFilter.mesh) < minDepth ? new List<Triangle>(triangleCandidateBoat).Find(p => p.TriangleIndex == currentTriangle.N2) : currentTriangle;
            currentTriangle = TriangleDepth(currentTriangle.N3, meshFilter.mesh) < minDepth ? new List<Triangle>(triangleCandidateBoat).Find(p => p.TriangleIndex == currentTriangle.N3) : currentTriangle;

            if (!new List<Triangle>(triangleCandidateBoat).Contains(currentTriangle))
            {
                stillInCandidates = false;
            }

            //We get all the triangle by graph traversal (BFS)
            Queue<Triangle> Q = new Queue<Triangle>();
            Q.Enqueue(currentTriangle);
            bottomHalf.Add(currentTriangle);

            while (Q.Count > 0)
            {
                currentTriangle = Q.Dequeue();
                foreach (int neighborIndex in currentTriangle.GetNeighborIndices())
                {
                    currentTriangle = new List<Triangle>(triangleCandidateBoat).Find(p => p.TriangleIndex == neighborIndex);
                    if (!new List<Triangle>(triangleCandidateBoat).Contains(currentTriangle) && !bottomHalf.Contains(currentTriangle))
                    {
                        bottomHalf.Add(currentTriangle);
                        Q.Enqueue(currentTriangle);
                    }
                }
            }
        }

        Triangle[] array = new Triangle[bottomHalf.Count];
        bottomHalf.CopyTo(array);

        return array;


    }

    private float TriangleDepth(int triangleIndex, Mesh mesh)
    {
        Vector3 v1 = mesh.vertices[mesh.triangles[triangleIndex * 3]];
        Vector3 v2 = mesh.vertices[mesh.triangles[triangleIndex * 3 + 1]];
        Vector3 v3 = mesh.vertices[mesh.triangles[triangleIndex * 3 + 2]];

        return ((v1 + v2 + v3) / 3f).y;
    }

    #endregion

    #region Compute forces on rigidbody ==============================================================================œ

    private void ApplyBuoyancy(Vector3[] verticesBoat, int[] trianglesBoat, Vector3[] verticesSea, int[] trianglesSea)
    {



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
                Gizmos.DrawSphere(transform.TransformPoint(chain[i + 1]), 0.01f);
                Gizmos.DrawLine(transform.TransformPoint(chain[i]), transform.TransformPoint(chain[i + 1]));
            }

        }
    }

    void DrawTriangleHelper(int i, Mesh mesh, Color color)
    {
        Vector3 v1 = mesh.vertices[mesh.triangles[i * 3]];
        Vector3 v2 = mesh.vertices[mesh.triangles[i * 3 + 1]];
        Vector3 v3 = mesh.vertices[mesh.triangles[i * 3 + 2]];

        Vector3 start = mesh == meshFilter.mesh ? Vector3.zero : Vector3.down * 2;
        Debug.DrawLine(start, (v1 + v2 + v3) / 3f, color, 5000f);
    }


    #endregion


}
