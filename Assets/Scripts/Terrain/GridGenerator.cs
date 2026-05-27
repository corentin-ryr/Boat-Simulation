using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    [RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
    public class GridGenerator : MonoBehaviour
    {
        public int seed;
        public bool randomizeSeed;

        [Header("Grid parameters")]
        public int nbIterRelaxation = 2;
        public bool normalizedRelaxation = false;
        public float hexagonRadius = 1f;
        public int gridSize = 10;

        [Header("Height map parameters")]
        public int numOctaves = 7;
        public float persistence = .5f;
        public float lacunarity = 2;
        public float initialScale = 2;
        public float elevationScale = 10f;

        [Header("Rendering")]
        public Material groundMaterial;

        [Header("Debug")]
        public bool showPrimalGizmos = false;
        public Color primalGizmoColor = new Color(0.3f, 0.8f, 1f);
        public Color edgeVertexColor = new Color(1f, 0.3f, 0.3f);
        public float edgeVertexSize = 0.08f;
        public bool showDualGizmos = false;
        public Color dualGizmoColor = new Color(1f, 0.8f, 0.2f);

        public HashSet<Polygon> DualPolygons { get; private set; } = new HashSet<Polygon>();
        public VertexCollection DualVertices { get; private set; } = new VertexCollection();

        List<Polygon> primalPolygons = new List<Polygon>();

        System.Random random;

        void Start()
        {
            seed = randomizeSeed ? UnityEngine.Random.Range(-10000, 10000) : seed;
            random = new System.Random(seed);

            VertexCollection primalVerts;
            (primalPolygons, primalVerts) = PolygonGridGenerator.GeneratePrimal(gridSize, hexagonRadius, random, Vector3.zero, nbIterRelaxation, normalizedRelaxation);
            (DualPolygons, DualVertices) = PolygonGridGenerator.GenerateDual(primalVerts);

            GenerateHeightMap();
            GenerateMesh();

            Debug.Log($"[GridGenerator] {DualPolygons.Count} polygons, {DualVertices.Count} vertices");
        }

        private void GenerateMesh()
        {
            List<Vector3> vertices = new List<Vector3>();
            List<int> triangles = new List<int>();

            foreach (Polygon polygon in DualPolygons)
            {
                Vector3[] verts = polygon.GetVerticesPosition();
                if (verts.Length < 3) continue;

                int baseIndex = vertices.Count;
                vertices.AddRange(verts);

                // Fan triangulation from vertex 0
                for (int i = 1; i < verts.Length - 1; i++)
                {
                    triangles.Add(baseIndex);
                    triangles.Add(baseIndex + i + 1);
                    triangles.Add(baseIndex + i);
                }
            }

            Mesh mesh = new Mesh();
            mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
            mesh.vertices = vertices.ToArray();
            mesh.triangles = triangles.ToArray();
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();

            GetComponent<MeshFilter>().mesh = mesh;
            if (groundMaterial != null)
                GetComponent<MeshRenderer>().material = groundMaterial;
        }

        private void GenerateHeightMap()
        {
            Vector2[] offsets = new Vector2[numOctaves];
            for (int i = 0; i < numOctaves; i++)
                offsets[i] = new Vector2(random.Next(-1000, 1000), random.Next(-1000, 1000));

            float minValue = float.MaxValue;
            float maxValue = float.MinValue;

            foreach (Vertex vertex in DualVertices.ToArray())
            {
                float x = vertex.Position.x;
                float z = vertex.Position.z;

                float noiseValue = 0;
                float scale = initialScale;
                float weight = 1;
                for (int i = 0; i < numOctaves; i++)
                {
                    Vector2 p = offsets[i] + new Vector2(x / gridSize, z / gridSize) * scale;
                    noiseValue += Mathf.PerlinNoise(p.x, p.y) * weight;
                    weight *= persistence;
                    scale *= lacunarity;
                }

                vertex.SetHeight(noiseValue);
                minValue = Mathf.Min(noiseValue, minValue);
                maxValue = Mathf.Max(noiseValue, maxValue);
            }

            if (maxValue != minValue)
            {
                foreach (Vertex vertex in DualVertices.ToArray())
                    vertex.SetHeight((vertex.Position.y - minValue) / (maxValue - minValue) * elevationScale);
            }
        }

        #region Debug

        void OnDrawGizmos()
        {
            if (showPrimalGizmos)
            {
                Gizmos.color = primalGizmoColor;
                foreach (Polygon p in primalPolygons)
                {
                    Vector3[] verts = p.GetVerticesPosition();
                    for (int i = 0; i < verts.Length; i++)
                        Gizmos.DrawLine(verts[i], verts[(i + 1) % verts.Length]);
                }

                Gizmos.color = edgeVertexColor;
                foreach (Polygon p in primalPolygons)
                    foreach (Vertex v in p.GetVertices())
                        if (v.IsEdge) Gizmos.DrawSphere(v.Position, edgeVertexSize);
            }

            if (showDualGizmos)
            {
                Gizmos.color = dualGizmoColor;
                foreach (Polygon p in DualPolygons)
                {
                    Vector3[] verts = p.GetVerticesPosition();
                    for (int i = 0; i < verts.Length; i++)
                        Gizmos.DrawLine(verts[i], verts[(i + 1) % verts.Length]);
                }
            }
        }

        public static void ShowMesh(Polygon[] polygons, Color color)
        {
            foreach (Polygon polygon in polygons)
            {
                Vector3[] verts = polygon.GetVerticesPosition();
                for (int i = 0; i < verts.Length; i++)
                    Debug.DrawLine(verts[i], verts[(i + 1) % verts.Length], color);
            }
        }

        public static void ShowVertices(Vertex[] vertices, Color color)
        {
            foreach (Vertex vertex in vertices)
            {
                foreach (Polygon polygon in vertex.Polygons)
                    Debug.DrawLine(vertex.Position, polygon.GetCenter(), color);
            }
        }

        #endregion
    }
}
