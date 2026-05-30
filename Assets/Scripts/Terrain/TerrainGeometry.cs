using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    // Coarse semantic classification of a primal cell, derived once at generation from the
    // elevation field. Lets gameplay code (NPCs, building placement, biome triggers) reason
    // about cells without re-sampling the noise field. Mirrors round-trip through ChunkSnapshot.
    //   Ocean   — every corner of the cell is at Y=0 (under the elevation threshold)
    //   Coastal — at least one corner at Y=0 and at least one above
    //   Land    — every corner above Y=0
    public enum CellTerrain : byte
    {
        Ocean = 0,
        Coastal = 1,
        Land = 2,
    }

    public class Polygon
    {
        Vertex[] vertices;
        Polygon[] neighbors;

        // Semantic class set by TerrainModel after vertex elevations are sampled. Default Ocean
        // is the correct value for a freshly-built flat (Y=0) chunk; classification overrides it.
        public CellTerrain Terrain = CellTerrain.Ocean;

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

        // Allocation-free adjacency setup, used by ComputeNeighbors' edge-map pass.
        public void InitNeighbors() { neighbors = new Polygon[vertices.Length]; }
        public void SetNeighborAt(int edge, Polygon neighbor) { neighbors[edge] = neighbor; }

        public Vector3[] GetVerticesPosition()
        {
            Vector3[] positions = new Vector3[vertices.Length];
            for (int i = 0; i < vertices.Length; i++) positions[i] = vertices[i].Position;
            return positions;
        }

        public Vertex[] GetVertices() => vertices;

        public Vertex[] GetEdge(int edgeNumber)
        {
            return new Vertex[] { GetVertices()[edgeNumber % vertices.Length], GetVertices()[(edgeNumber + 1) % vertices.Length] };
        }

        // Centroid, summed in place — no intermediate List/array/LINQ. This is called heavily
        // in relaxation (≈2× per polygon per iteration), so keeping it allocation-free matters.
        public Vector3 GetCenter()
        {
            Vector3 sum = Vector3.zero;
            for (int i = 0; i < vertices.Length; i++) sum += vertices[i].Position;
            return sum / vertices.Length;
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
        // The live incident-polygon set, exposed read-only so callers can enumerate it without
        // the per-access array allocation the old ToArray() getter caused. Do not mutate the
        // returned collection — use AddPolygon/RemovePolygon.
        public IReadOnlyCollection<Polygon> Polygons => polygons;
        public bool IsEdge { get => isEdge; set => isEdge = value; }

        public void AddPolygon(Polygon polygon) { polygons.Add(polygon); }
        public void RemovePolygon(Polygon polygon) { polygons.Remove(polygon); }

        public void SetHeight(float height)
        {
            position = new Vector3(position.x, height, position.z);
        }

        public void AccumulateMovement(Vector3 movement) { this.movement += movement; }

        // Unconditional move, used by the cross-chunk border relaxation pass which
        // deliberately moves IsEdge vertices (UpdateVertexPosition leaves them pinned).
        public void Translate(Vector3 delta) { position += delta; }

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

        // Round to 3 decimal places to absorb float drift across chunk offsets
        private static Vector3 Key(Vector3 p) => new Vector3(
            Mathf.Round(p.x * 1000f) / 1000f,
            Mathf.Round(p.y * 1000f) / 1000f,
            Mathf.Round(p.z * 1000f) / 1000f
        );

        public Vertex AddOrCreate(Vector3 position, bool isEdge = false)
        {
            Vector3 key = Key(position);
            if (key.x < minX) minX = key.x;
            if (key.x > maxX) maxX = key.x;
            if (key.z < minX) minY = key.z;
            if (key.z > maxX) maxY = key.z;

            if (!vertices.TryGetValue(key, out Vertex existing))
            {
                Vertex newVertex = new Vertex(position, isEdge);
                vertices[key] = newVertex;
                return newVertex;
            }
            return existing;
        }

        public Vertex GetAt(Vector3 position)
        {
            vertices.TryGetValue(Key(position), out Vertex v);
            return v;
        }

        public void AddVertex(Vertex vertex)
        {
            vertices[Key(vertex.Position)] = vertex;
        }

        public void Remove(Vertex vertex)
        {
            vertices.Remove(Key(vertex.Position));
        }

        public int Count { get => vertices.Count; }
        public Vertex[] ToArray() => vertices.Values.ToArray();

        // Allocation-free iteration over the contained vertices. Prefer this to ToArray()
        // anywhere the call site only enumerates and never indexes — saves the per-call
        // array allocation, which matters on the worker thread's BuildDual/DeepCopy loops.
        public Dictionary<Vector3, Vertex>.ValueCollection Values => vertices.Values;
        public float GetMapSize() => Mathf.Max(maxX - minX, maxY - minY);
    }
}
