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
                        // chain.Add(currentP);


                        goto End;
                    }
                }
            }
        }

        return (null, null);
    End:

        Vector3? vectorClosestToStart = null;
        //We have the first triangle, we start looping
        for (int i = 0; i < 100; i++) //Max number of iteration (security)
        {
            foreach (Triangle currentTriangle in nextCurrentTriangles) //Normally only one triangle in this array but if we have cases 3 or 4 we can have 0 or more than 1 triangles.
            {
                bool nextSwap;
                (currentCase, nextCurrentTriangles, nextSwap, nextTriangleToCheckAgainst) = IdentifyCase(currentTriangle, triangleToCheckAgainst, ref currentP, swaped);
                if (currentCase > 0)
                {
                    Triangle waterTriangle = swaped ? currentTriangle : triangleToCheckAgainst;
                    Triangle floatingTriangle = swaped ? triangleToCheckAgainst : currentTriangle;
                    bool changeFloatingTriangle = floatingTriangle != (nextSwap ? nextTriangleToCheckAgainst : nextCurrentTriangles[0]);

                    if (vectorClosestToStart is not null)
                    {
                        Vector3 waterNormal = Vector3.Cross(waterTriangle.Vertex1 - waterTriangle.Vertex2, waterTriangle.Vertex1 - waterTriangle.Vertex3).normalized;
                        Vector3 floatingNormal = Vector3.Cross(floatingTriangle.Vertex1 - floatingTriangle.Vertex2, floatingTriangle.Vertex1 - floatingTriangle.Vertex3).normalized;

                        List<Vector3> underwaterPoints = new List<Vector3>();

                        if (Vector3.Dot(floatingTriangle.Vertex1 - chain[chain.Count - 1], waterNormal) < 0) underwaterPoints.Add(floatingTriangle.Vertex1);
                        if (Vector3.Dot(floatingTriangle.Vertex2 - chain[chain.Count - 1], waterNormal) < 0) underwaterPoints.Add(floatingTriangle.Vertex2);
                        if (Vector3.Dot(floatingTriangle.Vertex3 - chain[chain.Count - 1], waterNormal) < 0) underwaterPoints.Add(floatingTriangle.Vertex3);

                        if (underwaterPoints.Count == 1 && FloaterHelper.showRays)
                        {
                            Vector3 tempNormal = Vector3.Cross(underwaterPoints[0] - chain[chain.Count - 1], underwaterPoints[0] - currentP);
                            if (Vector3.Dot(floatingNormal, tempNormal) > 0) lowerRing.Add(new Triangle(0, underwaterPoints[0], chain[chain.Count - 1], currentP));
                            else lowerRing.Add(new Triangle(0, chain[chain.Count - 1], underwaterPoints[0], currentP));

                        }
                        else if (underwaterPoints.Count == 2)
                        {
                            int indexStartingFloatingPoint = (vectorClosestToStart - underwaterPoints[0])?.magnitude < (vectorClosestToStart - underwaterPoints[1])?.magnitude ? 0 : 1;
                            bool rightSens = Vector3.Dot(floatingNormal, Vector3.Cross(underwaterPoints[indexStartingFloatingPoint] - chain[chain.Count - 1],
                                                                                        underwaterPoints[indexStartingFloatingPoint] - currentP)) > 0;

                            if (rightSens) lowerRing.Add(new Triangle(0, underwaterPoints[indexStartingFloatingPoint], chain[chain.Count - 1], currentP));
                            else lowerRing.Add(new Triangle(0, chain[chain.Count - 1], underwaterPoints[indexStartingFloatingPoint], currentP));

                            if (changeFloatingTriangle)
                            {
                                if (rightSens) lowerRing.Add(new Triangle(0, underwaterPoints[1 - indexStartingFloatingPoint], underwaterPoints[indexStartingFloatingPoint], currentP));
                                else lowerRing.Add(new Triangle(0, underwaterPoints[indexStartingFloatingPoint], underwaterPoints[1 - indexStartingFloatingPoint], currentP));
                            }
                        }
                        else if (FloaterHelper.showRays) Debug.Log("All three points underwater");

                    }

                    if (changeFloatingTriangle) // If we change the floating triangle at the next step
                    {
                        vectorClosestToStart = FloaterHelper.PointAdjacentOnEdge(nextSwap ? nextTriangleToCheckAgainst : nextCurrentTriangles[0], nextSwap ? nextCurrentTriangles[0] : nextTriangleToCheckAgainst, currentP);
                    }
                    if (vectorClosestToStart is not null) chain.Add(currentP);
                    swaped = nextSwap;
                    triangleToCheckAgainst = nextTriangleToCheckAgainst;
                    break; // No need to check the other potential triangles
                }
                UnityEngine.Debug.Log("No intersection");
                return (null, null);
            }


            if (chain.Count > 1 && (currentP - chain[0]).magnitude < 1E-3) break;

            if (i == 99) Debug.Log("Max number of iteration reached");
        }

        return (chain, lowerRing);
    }

    public enum PolygonState
    {
        start,
        middle,
        end
    }

    public static bool showRays = true;
    public static void HandlePolygons(Triangle triangle, Triangle waterTriangle, Vector3 pointA, Vector3 pointB, List<Triangle> polygons, Transform transform, bool newTriangle)
    {
        Debug.Log("Handle polygon");

        Vector3 waterNormal = Vector3.Cross(waterTriangle.Vertex1 - waterTriangle.Vertex2, waterTriangle.Vertex1 - waterTriangle.Vertex3).normalized;


        List<Vector3> underwaterPoints = new List<Vector3>();

        if (Vector3.Dot(triangle.Vertex1 - pointA, waterNormal) < 0) underwaterPoints.Add(triangle.Vertex1);
        if (Vector3.Dot(triangle.Vertex2 - pointA, waterNormal) < 0) underwaterPoints.Add(triangle.Vertex2);
        if (Vector3.Dot(triangle.Vertex3 - pointA, waterNormal) < 0) underwaterPoints.Add(triangle.Vertex3);

        if (underwaterPoints.Count == 1 && FloaterHelper.showRays)
        {
            polygons.Add(new Triangle(0, underwaterPoints[0], pointA, pointB));
            // Debug.DrawLine(transform.TransformPoint(underwaterPoints[0]), transform.TransformPoint(pointA), Color.green);
            // Debug.DrawLine(transform.TransformPoint(underwaterPoints[0]), transform.TransformPoint(pointB), Color.green);
            // Debug.DrawLine(transform.TransformPoint(pointA), transform.TransformPoint(pointB), Color.green);

            // FloaterHelper.showRays = false;
        }
        else if (underwaterPoints.Count == 2)
        {
            int indexStartingFloatingPoint = (pointA - underwaterPoints[0]).magnitude < (pointB - underwaterPoints[0]).magnitude ? 0 : 1;
            // Vector3 vectorClosestToStart = FloaterHelper.PointAdjacentOnEdge(triangle, waterTriangle, pointA);

            polygons.Add(new Triangle(0, underwaterPoints[indexStartingFloatingPoint], pointA, pointB));
            Debug.DrawLine(transform.TransformPoint(underwaterPoints[indexStartingFloatingPoint]), transform.TransformPoint(pointA), Color.green);
            Debug.DrawLine(transform.TransformPoint(underwaterPoints[indexStartingFloatingPoint]), transform.TransformPoint(pointB), Color.green);
            Debug.DrawLine(transform.TransformPoint(pointA), transform.TransformPoint(pointB));
            if (newTriangle)
            {
                Debug.DrawLine(transform.TransformPoint(underwaterPoints[1 - indexStartingFloatingPoint]), transform.TransformPoint(underwaterPoints[indexStartingFloatingPoint]), Color.red);
                Debug.DrawLine(transform.TransformPoint(underwaterPoints[1 - indexStartingFloatingPoint]), transform.TransformPoint(pointB), Color.red);
                Debug.DrawLine(transform.TransformPoint(underwaterPoints[indexStartingFloatingPoint]), transform.TransformPoint(pointB), Color.red);
                polygons.Add(new Triangle(0, underwaterPoints[1 - indexStartingFloatingPoint], underwaterPoints[indexStartingFloatingPoint], pointB));
            }
        }
        else if (FloaterHelper.showRays) Debug.Log("All three points underwater");

    }

    public static Vector3 PointAdjacentOnEdge(Triangle triangle, Triangle waterTriangle, Vector3 point)
    {
        float edge1 = Vector3.Cross((triangle.Vertex1 - triangle.Vertex2).normalized, triangle.Vertex1 - point).magnitude;
        float edge2 = Vector3.Cross((triangle.Vertex2 - triangle.Vertex3).normalized, (triangle.Vertex2 - point).normalized).magnitude;
        float edge3 = Vector3.Cross((triangle.Vertex3 - triangle.Vertex1).normalized, (triangle.Vertex3 - point).normalized).magnitude;
        float[] edges = new float[] { edge1, edge2, edge3 };

        int edgeIndex = edges.ToList().IndexOf(edges.Min());

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
