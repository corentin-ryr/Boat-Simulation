using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    public class ChunkManager : MonoBehaviour
    {
        [Header("Grid")]
        public int chunkGridSize = 5;
        public float hexRadius = 1f;
        public int nbIterRelaxation = 2;
        public bool normalizedRelaxation = false;
        public int seed = 0;

        [Header("Streaming")]
        public int loadRadius = 2;
        public Transform cameraTarget;

        [Header("Rendering")]
        public Material groundMaterial;

        readonly Dictionary<ChunkCoord, List<Polygon>> loadedPrimalChunks = new();
        readonly VertexCollection globalPrimalVerts = new();
        readonly Dictionary<ChunkCoord, GameObject> chunkObjects = new();

        ChunkCoord lastCameraChunk;

        ChunkCoord WorldToChunk(Vector3 pos) => ChunkCoord.FromWorldPos(pos, hexRadius, chunkGridSize);
        Vector3 ChunkToWorld(ChunkCoord coord) => coord.WorldCenter(hexRadius, chunkGridSize);

        void Start()
        {
            if (cameraTarget == null) cameraTarget = Camera.main?.transform;
            lastCameraChunk = WorldToChunk(cameraTarget.position);
            UpdateChunks(lastCameraChunk);
        }

        void Update()
        {
            ChunkCoord current = WorldToChunk(cameraTarget.position);
            if (current != lastCameraChunk)
            {
                lastCameraChunk = current;
                UpdateChunks(current);
            }
        }

        void UpdateChunks(ChunkCoord center)
        {
            HashSet<ChunkCoord> target = new(center.HexesInRange(loadRadius));
            HashSet<ChunkCoord> loaded = new(loadedPrimalChunks.Keys);

            foreach (ChunkCoord coord in loaded.Except(target).ToList())
                UnloadChunk(coord);

            foreach (ChunkCoord coord in target.Except(loaded))
                LoadChunk(coord);

            RebuildMeshes();
        }

        void LoadChunk(ChunkCoord coord)
        {
            int chunkSeed = seed ^ (coord.q * 73856093) ^ (coord.r * 19349663);
            var random = new System.Random(chunkSeed);

            var (polygons, verts) = PolygonGridGenerator.GeneratePrimal(
                chunkGridSize, hexRadius, random, ChunkToWorld(coord), nbIterRelaxation, normalizedRelaxation);

            PolygonGridGenerator.Stitch(polygons, verts, globalPrimalVerts);
            loadedPrimalChunks[coord] = polygons;
        }

        void UnloadChunk(ChunkCoord coord)
        {
            foreach (Polygon p in loadedPrimalChunks[coord])
            {
                foreach (Vertex v in p.GetVertices())
                {
                    v.RemovePolygon(p);
                    if (v.Polygons.Length == 0)
                        globalPrimalVerts.Remove(v);
                }
            }
            loadedPrimalChunks.Remove(coord);
        }

        void RebuildMeshes()
        {
            var (dualPolygons, _) = PolygonGridGenerator.GenerateDual(globalPrimalVerts);

            // Group dual polygons by which chunk their center falls in
            var byChunk = new Dictionary<ChunkCoord, List<Polygon>>();
            foreach (Polygon p in dualPolygons)
            {
                ChunkCoord coord = WorldToChunk(p.GetCenter());
                if (!byChunk.TryGetValue(coord, out var list))
                    byChunk[coord] = list = new List<Polygon>();
                list.Add(p);
            }

            // Build / update mesh objects
            foreach (var (coord, polygons) in byChunk)
            {
                if (!chunkObjects.TryGetValue(coord, out GameObject go))
                {
                    go = new GameObject($"Chunk {coord}");
                    go.transform.SetParent(transform, false);
                    go.AddComponent<MeshFilter>();
                    go.AddComponent<MeshRenderer>().material = groundMaterial;
                    chunkObjects[coord] = go;
                }
                BuildMesh(go.GetComponent<MeshFilter>(), polygons);
            }

            // Destroy chunk objects whose dual polygons are gone
            foreach (ChunkCoord coord in chunkObjects.Keys.Except(byChunk.Keys).ToList())
            {
                Destroy(chunkObjects[coord]);
                chunkObjects.Remove(coord);
            }
        }

        static void BuildMesh(MeshFilter mf, List<Polygon> polygons)
        {
            var vertices = new List<Vector3>();
            var triangles = new List<int>();

            foreach (Polygon p in polygons)
            {
                Vector3[] verts = p.GetVerticesPosition();
                if (verts.Length < 3) continue;
                int baseIndex = vertices.Count;
                vertices.AddRange(verts);
                for (int i = 1; i < verts.Length - 1; i++)
                {
                    triangles.Add(baseIndex);
                    triangles.Add(baseIndex + i + 1);
                    triangles.Add(baseIndex + i);
                }
            }

            Mesh mesh = new Mesh
            {
                indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
                vertices = vertices.ToArray(),
                triangles = triangles.ToArray()
            };
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();
            mf.mesh = mesh;
        }
    }
}
