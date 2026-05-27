using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    // The render layer: turns primal chunks (from a TerrainModel) into dual grids and Unity
    // meshes. A pure consumer — it reads the model and never writes to it. Meshes are built
    // only for the coords it is asked to show, and freed when no longer wanted, so the
    // expensive dual/mesh work tracks the camera while the primal data persists.
    public class ChunkRenderer
    {
        class RenderedChunk
        {
            public GameObject Go;
            public List<Polygon> Dual;
            public int BuiltFromVersion;
        }

        readonly TerrainModel model;
        readonly Transform parent;
        readonly Material material;

        readonly Dictionary<ChunkCoord, RenderedChunk> rendered = new();

        public ChunkRenderer(TerrainModel model, Transform parent, Material material)
        {
            this.model = model;
            this.parent = parent;
            this.material = material;
        }

        // Dual polygons of every currently rendered chunk (for debug gizmos).
        public IEnumerable<List<Polygon>> DualPolygons => rendered.Values.Select(r => r.Dual);

        // Render exactly the chunks in `desired`: build/refresh those, free the rest.
        public void SyncTo(IEnumerable<ChunkCoord> desired)
        {
            HashSet<ChunkCoord> want = new HashSet<ChunkCoord>(desired);

            foreach (ChunkCoord coord in rendered.Keys.Where(c => !want.Contains(c)).ToList())
                Hide(coord);

            foreach (ChunkCoord coord in want)
                Show(coord);
        }

        void Show(ChunkCoord coord)
        {
            // Skip if the primal isn't in the mirror yet (it streams in asynchronously; the
            // render set's halo means neighbours are normally present, but stay safe).
            if (!model.TryGet(coord, out PrimalChunk primal)) return;

            if (rendered.TryGetValue(coord, out RenderedChunk rc))
            {
                if (rc.BuiltFromVersion == primal.Version) return; // mesh already up to date
            }
            else
            {
                GameObject go = new GameObject($"Chunk {coord}");
                go.transform.SetParent(parent, false);
                go.AddComponent<MeshFilter>();
                go.AddComponent<MeshRenderer>().material = material;
                rc = new RenderedChunk { Go = go };
                rendered[coord] = rc;
            }

            BuildNeighborBorderMaps(coord, out var neighborFaces, out var minNeighborCoord);

            var (dualPolygons, _) = PolygonGridGenerator.GenerateDual(primal, model.hexRadius, neighborFaces, minNeighborCoord);
            rc.Dual = dualPolygons.ToList();
            rc.BuiltFromVersion = primal.Version;

            BuildMesh(rc.Go.GetComponent<MeshFilter>(), rc.Dual);
        }

        // Gather, from the chunk's 6 loaded neighbours, the faces incident to each shared
        // border vertex (for cell completion) and the lowest neighbour coord per border
        // position (for ownership). Keyed by PolygonGridGenerator.LatticeKey.
        void BuildNeighborBorderMaps(ChunkCoord coord,
            out Dictionary<(int, int), List<Polygon>> neighborFaces,
            out Dictionary<(int, int), ChunkCoord> minNeighborCoord)
        {
            neighborFaces = new Dictionary<(int, int), List<Polygon>>();
            minNeighborCoord = new Dictionary<(int, int), ChunkCoord>();

            foreach (ChunkCoord n in coord.HexesInRange(1))
            {
                if (n == coord) continue;
                if (!model.TryGet(n, out PrimalChunk np)) continue;

                foreach (Vertex v in np.Verts.ToArray())
                {
                    if (!v.IsEdge) continue;
                    (int, int) key = PolygonGridGenerator.LatticeKey(v.Position, model.hexRadius);

                    if (!neighborFaces.TryGetValue(key, out List<Polygon> list))
                        neighborFaces[key] = list = new List<Polygon>();
                    list.AddRange(v.Polygons);

                    if (!minNeighborCoord.TryGetValue(key, out ChunkCoord cur) || n.CompareTo(cur) < 0)
                        minNeighborCoord[key] = n;
                }
            }
        }

        void Hide(ChunkCoord coord)
        {
            if (!rendered.TryGetValue(coord, out RenderedChunk rc)) return;
            Object.Destroy(rc.Go);
            rendered.Remove(coord);
        }

        public void Clear()
        {
            foreach (RenderedChunk rc in rendered.Values) Object.Destroy(rc.Go);
            rendered.Clear();
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
