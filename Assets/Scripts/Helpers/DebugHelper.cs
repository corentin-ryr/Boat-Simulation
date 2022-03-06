using System.Collections;
using System.Collections.Generic;
using Habrador_Computational_Geometry;
using UnityEngine;

public static class DebugHelper
{

    public static void DrawBounds(Bounds b, float delay = 0, GameObject objectReference = null)
    {
        Vector3 offset = objectReference?.transform.position ?? Vector3.zero;
        // bottom
        var p1 = new Vector3(b.min.x, b.min.y, b.min.z) + offset;
        var p2 = new Vector3(b.max.x, b.min.y, b.min.z) + offset;
        var p3 = new Vector3(b.max.x, b.min.y, b.max.z) + offset;
        var p4 = new Vector3(b.min.x, b.min.y, b.max.z) + offset;

        Debug.DrawLine(p1, p2, Color.blue, delay);
        Debug.DrawLine(p2, p3, Color.red, delay);
        Debug.DrawLine(p3, p4, Color.yellow, delay);
        Debug.DrawLine(p4, p1, Color.magenta, delay);

        // top
        var p5 = new Vector3(b.min.x, b.max.y, b.min.z) + offset;
        var p6 = new Vector3(b.max.x, b.max.y, b.min.z) + offset;
        var p7 = new Vector3(b.max.x, b.max.y, b.max.z) + offset;
        var p8 = new Vector3(b.min.x, b.max.y, b.max.z) + offset;

        Debug.DrawLine(p5, p6, Color.blue, delay);
        Debug.DrawLine(p6, p7, Color.red, delay);
        Debug.DrawLine(p7, p8, Color.yellow, delay);
        Debug.DrawLine(p8, p5, Color.magenta, delay);

        // sides
        Debug.DrawLine(p1, p5, Color.white, delay);
        Debug.DrawLine(p2, p6, Color.gray, delay);
        Debug.DrawLine(p3, p7, Color.green, delay);
        Debug.DrawLine(p4, p8, Color.cyan, delay);
    }

    public static void ShowMesh(Triangle[] triangles, Transform transform, Color color)
    {
        Debug.Log("Number of candidates: " + triangles.Length);
        Vector3[] debugVertices = new Vector3[triangles.Length * 3];
        Matrix4x4 transformation = transform.localToWorldMatrix;
        for (int i = 0; i < triangles.Length; i++)
        {
            debugVertices[i * 3] = transformation.MultiplyPoint3x4(triangles[i].Vertex1);
            debugVertices[i * 3 + 1] = transformation.MultiplyPoint3x4(triangles[i].Vertex2);
            debugVertices[i * 3 + 2] = transformation.MultiplyPoint3x4(triangles[i].Vertex3);

            Debug.DrawLine(debugVertices[i * 3], debugVertices[i * 3 + 1], color);
            Debug.DrawLine(debugVertices[i * 3 + 1], debugVertices[i * 3 + 2], color);
            Debug.DrawLine(debugVertices[i * 3 + 2], debugVertices[i * 3], color);
        }
    }
    public static void ShowMesh(Triangle[] triangles, Color color)
    {
        for (int i = 0; i < triangles.Length; i++)
        {
            Vector3 a = triangles[i].Vertex1;
            Vector3 b = triangles[i].Vertex2;
            Vector3 c = triangles[i].Vertex3;

            Debug.DrawLine(a, b, color);
            Debug.DrawLine(b, c, color);
            Debug.DrawLine(c, a, color);
        }
    }

    public static void ShowMesh(Vector3[] vertices, int[] triangles, Transform parentTransfrom, Color color)
    {
        for (int i = 0; i < triangles.Length / 3; i++)
        {
            Vector3 a = parentTransfrom.TransformPoint(vertices[triangles[i * 3]]);
            Vector3 b = parentTransfrom.TransformPoint(vertices[triangles[i * 3 + 1]]);
            Vector3 c = parentTransfrom.TransformPoint(vertices[triangles[i * 3 + 2]]);

            Debug.DrawLine(a, b, color);
            Debug.DrawLine(b, c, color);
            Debug.DrawLine(c, a, color);
        }
    }
 

}
