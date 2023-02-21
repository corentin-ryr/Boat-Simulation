using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public static class FloaterHelper
{
    //Triangle1 is from the boat and triangle2 is from the water
    //TODO compute cases 3 and 4
    public static (int, Triangle[], bool, Triangle) IdentifyCase(Triangle triangle1, Triangle triangle2, ref Vector3 P, bool swaped, Transform debugTransform = null)
    {
        Vector3 a = triangle1.Vertex1;
        Vector3 b = triangle1.Vertex2;
        Vector3 c = triangle1.Vertex3;
        Vector3 d = triangle2.Vertex1;
        Vector3 e = triangle2.Vertex2;
        Vector3 f = triangle2.Vertex3;

        //Case 1 if p is inside of triangle1 
        Triangle[] triangleIndicesToCheck1 = null;
        Vector3? case1Vec;
        if ((case1Vec = IntersectionEdgeTriangle(a, b, c, d, e, P)) is not null)
        {
            triangleIndicesToCheck1 = new Triangle[] { triangle2.T1 };
        }
        else if ((case1Vec = IntersectionEdgeTriangle(a, b, c, e, f, P)) is not null)
        {
            triangleIndicesToCheck1 = new Triangle[] { triangle2.T2 };
        }
        else if ((case1Vec = IntersectionEdgeTriangle(a, b, c, f, d, P)) is not null)
        {
            triangleIndicesToCheck1 = new Triangle[] { triangle2.T3 };
        }


        //Case 2 if p is inside of triangle2
        Triangle[] triangleIndicesToCheck2 = null;
        Vector3? case2Vec = null;
        if ((case2Vec = IntersectionEdgeTriangle(d, e, f, a, b, P)) is not null)
        {
            triangleIndicesToCheck2 = new Triangle[] { triangle1.T1 };
        }
        else if ((case2Vec = IntersectionEdgeTriangle(d, e, f, b, c, P)) is not null)
        {
            triangleIndicesToCheck2 = new Triangle[] { triangle1.T2 };
        }
        else if ((case2Vec = IntersectionEdgeTriangle(d, e, f, c, a, P)) is not null)
        {
            triangleIndicesToCheck2 = new Triangle[] { triangle1.T3 };
        }

        Triangle nextTriangleToCheckAgainst;

        //Case 3 if p lies in the other triangle
        if (case2Vec is not null)
        {
            if (case2Vec == a || case2Vec == b || case2Vec == c)
            {
                P = (Vector3)case2Vec;
                triangleIndicesToCheck1 = FindTriangleCommonPoint(triangle1, P);
                nextTriangleToCheckAgainst = triangle2;
                return (3, triangleIndicesToCheck1, swaped, nextTriangleToCheckAgainst);
            }
        }

        if (case1Vec is not null && case2Vec is not null && case1Vec == case2Vec) //Case when we have the point on the edge of the two triangles.
        {
            nextTriangleToCheckAgainst = triangleIndicesToCheck2[0];
            swaped = !swaped;
            P = (Vector3)case1Vec;
            return (5, triangleIndicesToCheck1, swaped, nextTriangleToCheckAgainst); //We act like a case 1 but we have a special triangleToCheck
        }

        if (case1Vec is not null) //Case 1, we have to swap and we set the triangle to check next
        {
            swaped = !swaped;
            nextTriangleToCheckAgainst = triangle1;
            P = (Vector3)case1Vec;
            return (1, triangleIndicesToCheck1, swaped, nextTriangleToCheckAgainst);
        }

        if (case2Vec is not null)//Case 2, we don't swap and we also set the next triangle to check
        {
            P = (Vector3)case2Vec;
            nextTriangleToCheckAgainst = triangle2;
            return (2, triangleIndicesToCheck2, swaped, nextTriangleToCheckAgainst);
        }

        return (0, triangleIndicesToCheck1, swaped, null);
    }

    public static Vector3? IntersectionEdgeTriangle(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 e, Vector3 previousP)
    {
        Vector3 n = Vector3.Cross(a - b, a - c).normalized;

        float E = Vector3.Dot(a - e, n);
        float D = Vector3.Dot(a - d, n);
        //Are d and e on different side of triangle abc ?
        if (D * E > 0) return null;

        float t = D / (D - E);
        Vector3 p = t * e + (1 - t) * d;
        if (p == previousP) return null;

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

    // For the case 3 we need to find all the triangles that touch the point "commonPoint".
    private static Triangle[] FindTriangleCommonPoint(Triangle startTriangle, Vector3 commonPoint)
    {
        List<Triangle> touchingTriangles = new List<Triangle>();

        //Find a point underwater
        Queue<Triangle> Q = new Queue<Triangle>();
        List<Triangle> triangleVisited = new List<Triangle>();
        Q.Enqueue(startTriangle);
        touchingTriangles.Add(startTriangle);

        while (Q.Count > 0)
        {
            Triangle currentTriangle = Q.Dequeue();
            foreach (Triangle neighbor in currentTriangle.GetNeighbors())
            {
                if (triangleVisited.Contains(neighbor)) continue;

                triangleVisited.Add(neighbor);

                if (neighbor.Vertex1 == commonPoint || neighbor.Vertex2 == commonPoint || neighbor.Vertex3 == commonPoint)
                {
                    touchingTriangles.Add(neighbor);
                    Q.Enqueue(neighbor);
                }

            }
        }

        return touchingTriangles.ToArray();
    }

    private static (Triangle[], bool, Triangle, Vector3) FindInitialIntersection(Cell[] gridCells)
    {
        Vector3 currentP = Vector3.positiveInfinity;
        bool swaped = false;
        Triangle triangleToCheckAgainst = null;
        Triangle[] nextCurrentTriangles = new Triangle[0];
        int currentCase;

        foreach (Cell cell in gridCells)
        {
            if (!cell.HasCandidates()) continue;
            foreach (Triangle triangleSet2 in cell.TriangleSet2)
            {
                foreach (Triangle triangleSet1 in cell.TriangleSet1)
                {
                    (currentCase, nextCurrentTriangles, swaped, triangleToCheckAgainst) = IdentifyCase(triangleSet1, triangleSet2, ref currentP, swaped);
                    if (currentCase > 0)
                    {
                        return (nextCurrentTriangles, swaped, triangleToCheckAgainst, currentP);
                    }
                }
            }
        }

        throw new System.Exception("No intersection");
    }



    public static (List<Vector3>, List<Triangle>) ComputeIntersectingLine(Cell[] gridCells, Transform debugTransform)
    {
        List<Vector3> chain = new List<Vector3>();

        //Find first intersection
        bool swaped = false;
        Triangle triangleToCheckAgainst = null;
        Vector3 currentP = Vector3.positiveInfinity;
        Triangle[] nextCurrentTriangles = new Triangle[0];
        Triangle[] currentTriangles = new Triangle[0];
        Triangle nextTriangleToCheckAgainst = null;

        List<Triangle> lowerRing = new List<Triangle>();

        try
        {
            (nextCurrentTriangles, swaped, triangleToCheckAgainst, currentP) = FloaterHelper.FindInitialIntersection(gridCells);
        }
        catch (System.Exception)
        {
            Debug.Log("No intersection");
            return (null, null);
        }

        currentTriangles = nextCurrentTriangles;


        Vector3? vectorClosestToStart = null;
        //We have the first triangle, we start looping
        for (int i = 0; i < 1000; i++) //Max number of iteration (security)
        {
            Color color = Color.Lerp(Color.black, Color.white, (float)i / (float)(chain.Count - 1));
            foreach (Triangle currentTriangle in currentTriangles)
            {
                bool nextSwap;
                int currentCase;

                (currentCase, nextCurrentTriangles, nextSwap, nextTriangleToCheckAgainst) = IdentifyCase(currentTriangle, triangleToCheckAgainst, ref currentP, swaped, debugTransform);
                if (currentCase > 0)
                {
                    Triangle waterTriangle = swaped ? currentTriangle : triangleToCheckAgainst;
                    Triangle floatingTriangle = swaped ? triangleToCheckAgainst : currentTriangle;
                    bool changeFloatingTriangle = floatingTriangle != (nextSwap ? nextTriangleToCheckAgainst : nextCurrentTriangles[0]);

                    // Compute the lower ring
                    if (vectorClosestToStart is not null)
                    {
                        Vector3 waterNormal = Vector3.Cross(waterTriangle.Vertex1 - waterTriangle.Vertex2, waterTriangle.Vertex1 - waterTriangle.Vertex3).normalized;
                        Vector3 floatingNormal = Vector3.Cross(floatingTriangle.Vertex1 - floatingTriangle.Vertex2, floatingTriangle.Vertex1 - floatingTriangle.Vertex3).normalized;

                        List<Vector3> underwaterPoints = new List<Vector3>();

                        if (Vector3.Dot(floatingTriangle.Vertex1 - chain[chain.Count - 1], waterNormal) < 0) underwaterPoints.Add(floatingTriangle.Vertex1);
                        if (Vector3.Dot(floatingTriangle.Vertex2 - chain[chain.Count - 1], waterNormal) < 0) underwaterPoints.Add(floatingTriangle.Vertex2);
                        if (Vector3.Dot(floatingTriangle.Vertex3 - chain[chain.Count - 1], waterNormal) < 0) underwaterPoints.Add(floatingTriangle.Vertex3);

                        if (underwaterPoints.Count == 1)
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
                        else Debug.Log("All three points underwater");

                    }
                    if (changeFloatingTriangle)
                    {
                        vectorClosestToStart = FloaterHelper.PointAdjacentOnEdge(nextSwap ? nextTriangleToCheckAgainst : nextCurrentTriangles[0], nextSwap ? nextCurrentTriangles[0] : nextTriangleToCheckAgainst, currentP);
                    }
                    if (vectorClosestToStart is not null) chain.Add(currentP);

                    swaped = nextSwap;
                    triangleToCheckAgainst = nextTriangleToCheckAgainst;
                    currentTriangles = nextCurrentTriangles;

                    break; // No need to check the other potential triangles
                }
            }

            if (nextCurrentTriangles is null) throw new System.Exception("Error during the intersection computation");
            if (chain.Count > 1 && (currentP == chain[0])) break;
            if (i == 999) Debug.Log("Max number of iteration reached");
        }

        return (chain, lowerRing);
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
