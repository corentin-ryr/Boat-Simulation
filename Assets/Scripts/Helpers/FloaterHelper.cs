using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public static class FloaterHelper
{
    //Triangle1 is from the boat and triangle2 is from the water
    //TODO compute cases 3 and 4
    public static (int, Triangle[], bool, Triangle) IdentifyCase(Triangle triangle1, Triangle triangle2, ref Vector3 P, bool swaped)
    {
        Vector3 a = triangle1.Vertex1;
        Vector3 b = triangle1.Vertex2;
        Vector3 c = triangle1.Vertex3;
        Vector3 d = triangle2.Vertex1;
        Vector3 e = triangle2.Vertex2;
        Vector3 f = triangle2.Vertex3;

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

        Triangle nextTriangleToCheckAgainst;

        if (case1 && case2 && tempP1 == tempP2) //Case when we have the point on the edge of the two triangles. The next triangle to check 
        {
            // triangleToCheck = swaped ? triangleIndicesToCheck2[0] : triangleIndicesToCheck1[0];
            nextTriangleToCheckAgainst = triangleIndicesToCheck2[0];
            swaped = !swaped;
            P = tempP1;
            return (1, triangleIndicesToCheck1, swaped, nextTriangleToCheckAgainst); //We act like a case 1 but we have a special triangleToCheck
        }
        else if (case1) //Case 1, we have to swap and we set the triangle to check next
        {
            swaped = !swaped;
            nextTriangleToCheckAgainst = triangle1;
            P = tempP1;
            return (1, triangleIndicesToCheck1, swaped, nextTriangleToCheckAgainst);
        }
        else if (case2)//Case 2, we don't swap and we also set the next triangle to check
        {
            P = tempP2;
            nextTriangleToCheckAgainst = triangle2;
            return (2, triangleIndicesToCheck2, swaped, nextTriangleToCheckAgainst);
        }

        return (0, triangleIndicesToCheck1, swaped, null);
    }

    public static bool IntersectionEdgeTriangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 e, ref Vector3 newP, Vector3 previousP)
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



    public static (List<Vector3>, List<Triangle>) ComputeIntersectingLine(Cell[] gridCells, Transform transform)
    {
        List<Vector3> chain = new List<Vector3>();

        //Find first intersection
        bool swaped = false;
        Triangle currentTriangle;
        Triangle triangleToCheckAgainst = null;
        Vector3 currentP = Vector3.positiveInfinity;
        Triangle[] nextCurrentTriangles = new Triangle[0];
        Triangle nextTriangleToCheckAgainst = null;
        int currentCase;

        List<Triangle> lowerRing = new List<Triangle>();

        foreach (Cell cell in gridCells)
        {
            foreach (Triangle triangleSet2 in cell.TriangleSet2)
            {
                foreach (Triangle triangleSet1 in cell.TriangleSet1)
                {
                    (currentCase, nextCurrentTriangles, swaped, triangleToCheckAgainst) = IdentifyCase(triangleSet1, triangleSet2, ref currentP, swaped);
                    if (currentCase > 0)
                    {
                        chain.Add(currentP);

                        goto End;
                    }
                }
            }
        }

        return (null, null);
    End:

        //We have the first triangle, we start looping
        for (int i = 0; i < 100; i++) //Max number of iteration (security)
        {
            foreach (Triangle nextCurrentTriangle in nextCurrentTriangles) //Normally only one triangle in this array but if we have cases 3 or 4 we can have 0 or more than 1 triangles.
            {
                currentTriangle = nextCurrentTriangle;
                bool nextSwap;
                (currentCase, nextCurrentTriangles, nextSwap, nextTriangleToCheckAgainst) = IdentifyCase(currentTriangle, triangleToCheckAgainst, ref currentP, swaped);
                if (currentCase > 0)
                {

                    // HandlePolygons(swaped ? triangleToCheckAgainst : currentTriangle, swaped ? currentTriangle : triangleToCheckAgainst, chain[chain.Count - 1], currentP, lowerRing, transform);

                    chain.Add(currentP);
                    swaped = nextSwap;
                    triangleToCheckAgainst = nextTriangleToCheckAgainst;
                    break;
                }
                UnityEngine.Debug.Log("No intersection");
                return (null, null);
            }


            if ((currentP - chain[0]).magnitude < 1E-3) break;

            if (i == 99) Debug.Log("Max number of iteration reached");
        }

        return (chain, null);
    }

    public enum PolygonState
    {
        start,
        middle,
        end
    }
    public static void HandlePolygons(Triangle triangle, Triangle waterTriangle, Vector3 pointA, Vector3 pointB, List<Triangle> polygons, Transform transform)
    {
        Debug.Log("Handle polygon");

        Vector3 waterNormal = Vector3.Cross(waterTriangle.Vertex1 - waterTriangle.Vertex2, waterTriangle.Vertex1 - waterTriangle.Vertex3).normalized;

        Debug.DrawRay(transform.TransformPoint(waterTriangle.Vertex1), transform.TransformDirection(waterNormal), Color.red);

        List<Vector3> underwaterPoints = new List<Vector3>();
        if (Vector3.Dot(triangle.Vertex1 - pointA, waterNormal) < 0) underwaterPoints.Add(triangle.Vertex1);
        if (Vector3.Dot(triangle.Vertex2 - pointA, waterNormal) < 0) underwaterPoints.Add(triangle.Vertex2);
        if (Vector3.Dot(triangle.Vertex3 - pointA, waterNormal) < 0) underwaterPoints.Add(triangle.Vertex3);
        
        if (underwaterPoints.Count == 1) {
            polygons.Add(new Triangle(0, underwaterPoints[0], pointA, pointB));
        }
        if (underwaterPoints.Count == 2) {
            // polygons.Add(new Triangle(0, underwaterPoints[1], pointA, pointB));
            // polygons.Add(new Triangle(0, underwaterPoints[0], pointA, pointB));
        }
        else Debug.Log("All three points underwater");
        
        // Get points in triangle that are under water

        // if (polygons.ContainsKey(triangle))
        // {
        //     if (state == PolygonState.start)
        //     {
        //         //Get the point underwater on the intersected edge
        //         polygons[triangle].Add(PointAdjacentOnEdge(triangle, waterTriangle, newPoint));
        //     }
        //     polygons[triangle].Add(newPoint);
        //     if (state == PolygonState.end)
        //     {
        //         //Get the point underwater on the intersected edge
        //         polygons[triangle].Add(PointAdjacentOnEdge(triangle, waterTriangle, newPoint));
        //     }
        // }
        // else if (state == PolygonState.start)
        // {
        //     List<Vector3> tempList = new List<Vector3>();
        //     tempList.Add(PointAdjacentOnEdge(triangle, waterTriangle, newPoint));
        //     tempList.Add(newPoint);
        //     polygons.Add(triangle, tempList);
        // }
    }

    public static Vector3 PointAdjacentOnEdge(Triangle triangle, Triangle waterTriangle, Vector3 point)
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

}