using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    public class Polygon
    {
        Vertex[] vertices;
        Polygon[] neighbors;

        public Polygon(Vertex[] _vertices)
        {
            vertices = _vertices;
            foreach (Vertex vertex in _vertices)
                vertex.AddPolygon(this);
        }

        public void SetVertex(Vertex vertex, int vertexIndex)
        {
            vertices[vertexIndex] = vertex;
            vertex.AddPolygon(this);
        }

        public void SetNeighbors(Polygon[] _neighbors) { neighbors = _neighbors; }

        public Vector3[] GetVerticesPosition()
        {
            List<Vector3> verticesPosition = new List<Vector3>();
            foreach (Vertex vertex in vertices)
                verticesPosition.Add(vertex.Position);
            return verticesPosition.ToArray();
        }

        public Vertex[] GetVertices() => vertices;

        public Vertex[] GetEdge(int edgeNumber)
        {
            return new Vertex[] { GetVertices()[edgeNumber % vertices.Length], GetVertices()[(edgeNumber + 1) % vertices.Length] };
        }

        public Vector3 GetCenter()
        {
            return GetVerticesPosition().Aggregate(new Vector3(0, 0, 0), (s, v) => s + v) / GetVerticesPosition().Length;
        }

        public Polygon[] GetNeighbors() => neighbors;

        public Polygon AdjacentPolygon(List<Polygon> otherPolygons)
        {
            foreach (Polygon polygon in otherPolygons)
            {
                HashSet<Vertex> thisVertices = new HashSet<Vertex>(GetVertices());
                HashSet<Vertex> otherVertices = new HashSet<Vertex>(polygon.GetVertices());
                if (thisVertices.Intersect(otherVertices).Count() == 2) return polygon;
            }
            return null;
        }

        public Polygon NeighborAtEdge(int edgeNumber)
        {
            Vertex[] edgeVertices = GetEdge(edgeNumber);
            foreach (Polygon polygon in GetNeighbors())
            {
                HashSet<Vertex> otherVertices = new HashSet<Vertex>(polygon.GetVertices());
                if (otherVertices.Contains(edgeVertices[0]) && otherVertices.Contains(edgeVertices[1])) return polygon;
            }
            return null;
        }

        public void UpdateNeigbhor(Polygon previousNeighbor, Polygon newNeighbor)
        {
            for (int i = 0; i < neighbors.Length; i++)
            {
                if (neighbors[i] == previousNeighbor) { neighbors[i] = newNeighbor; return; }
            }
        }
    }


    public class Vertex
    {
        Vector3 position;
        HashSet<Polygon> polygons = new HashSet<Polygon>();
        bool isEdge;
        Vector3 movement = Vector3.zero;

        public Vertex(Vector3 _position, bool _isEdge = false) { position = _position; isEdge = _isEdge; }

        public Vector3 Position { get => position; }
        public Polygon[] Polygons { get => polygons.ToArray(); }
        public bool IsEdge { get => isEdge; set => isEdge = value; }

        public void AddPolygon(Polygon polygon) { polygons.Add(polygon); }
        public void RemovePolygon(Polygon polygon) { polygons.Remove(polygon); }

        public void SetHeight(float height)
        {
            position = new Vector3(position.x, height, position.z);
        }

        public void AccumulateMovement(Vector3 movement) { this.movement += movement; }

        public void UpdateVertexPosition()
        {
            if (!isEdge) position += movement * 0.1f;
            movement = Vector3.zero;
        }
    }


    public class VertexCollection
    {
        Dictionary<Vector3, Vertex> vertices;

        private float maxX = float.MinValue;
        private float minX = float.MinValue;
        private float maxY = float.MinValue;
        private float minY = float.MinValue;

        public VertexCollection() { vertices = new Dictionary<Vector3, Vertex>(); }

        public Vertex AddOrCreate(Vector3 position, bool isEdge = false)
        {
            if (position.x < minX) minX = position.x;
            if (position.x > maxX) maxX = position.x;
            if (position.z < minX) minY = position.z;
            if (position.z > maxX) maxY = position.z;

            if (!vertices.TryGetValue(position, out Vertex existing))
            {
                Vertex newVertex = new Vertex(position, isEdge);
                vertices[position] = newVertex;
                return newVertex;
            }
            return existing;
        }

        public int Count { get => vertices.Count; }
        public Vertex[] ToArray() => vertices.Values.ToArray();
        public float GetMapSize() => Mathf.Max(maxX - minX, maxY - minY);
    }
}
