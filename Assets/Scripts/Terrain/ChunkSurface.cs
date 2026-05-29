using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    // Owns the main-thread "presence" of terrain chunks in the scene: one GameObject per
    // chunk any client wants, with a MeshRenderer when someone wants it drawn and a
    // MeshCollider when someone wants it physical.
    //
    // Multi-client by design: ChunkManager registers a render request (camera-rendered set)
    // *and* a collider request (so the boat doesn't fall through visible terrain). The
    // SimulationManager registers a collider-only request for the chunks under its agents.
    // Each client is expected to have already expanded its user-facing set to include the
    // minimal halo of lower-coord neighbours that own the perimeter cells — see
    // ChunkCoord.MinimalHalo. The surface itself takes the set at face value.
    //
    // Dual ownership is deterministic: each border cell is owned by the lowest-coord
    // contributing chunk, independent of which chunks happen to be presented. The dual is
    // cached on the PrimalChunk (so it survives across Apply calls, presentation toggles,
    // and per-frame churn). It's invalidated when the chunk's version changes (handled by
    // PrimalChunk.Version itself) or any neighbour's version changes (handled by the
    // cascade in TerrainModel.Install). The Unity Mesh is per-presence (Unity objects can't
    // travel) but rebuilt only when the cached dual has changed.
    public class ChunkSurface
    {
        class Presence
        {
            public GameObject Go;
            public MeshFilter Mf;
            public MeshRenderer Mr;     // null when not currently rendered
            public MeshCollider Mc;     // null when not currently collided
            public Mesh Mesh;
            public List<Polygon> DualRef; // last dual list we built the Mesh from (identity check)
            public int MeshBuiltFromVersion = -1;
        }

        readonly Transform parent;
        readonly Material material;
        readonly float hexRadius;

        readonly List<Func<ChunkCoord, PrimalChunk>> primalSources = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> renderRequests = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> colliderRequests = new();
        readonly Dictionary<ChunkCoord, Presence> presences = new();

        public ChunkSurface(Transform parent, Material material, float hexRadius)
        {
            this.parent = parent;
            this.material = material;
            this.hexRadius = hexRadius;
        }

        // For debug gizmos: the cached duals of every presented chunk's primal.
        public IEnumerable<List<Polygon>> DualPolygons
        {
            get
            {
                foreach (ChunkCoord coord in presences.Keys)
                {
                    PrimalChunk primal = LookupPrimal(coord);
                    if (primal?.Dual != null) yield return primal.Dual;
                }
            }
        }

        // Anyone with primal chunks the surface might need can register a lookup. Sources
        // are tried in registration order; first non-null wins. ChunkManager registers its
        // render mirror; SimulationManager registers its sim mirror as a fallback.
        public void AddPrimalSource(Func<ChunkCoord, PrimalChunk> source)
        {
            if (source != null) primalSources.Add(source);
        }

        // Replace this client's render set. Pass null to clear.
        public void SetRenderRequest(object client, IEnumerable<ChunkCoord> coords)
        {
            if (coords == null) renderRequests.Remove(client);
            else renderRequests[client] = new HashSet<ChunkCoord>(coords);
        }

        // Replace this client's collider set. Pass null to clear.
        public void SetColliderRequest(object client, IEnumerable<ChunkCoord> coords)
        {
            if (coords == null) colliderRequests.Remove(client);
            else colliderRequests[client] = new HashSet<ChunkCoord>(coords);
        }

        // Reconcile scene presence with the current requests. Cheap when nothing changed;
        // safe to call every frame. Steady-state cost is one HashSet build + iterate.
        public void Apply()
        {
            HashSet<ChunkCoord> renderUnion = Union(renderRequests);
            HashSet<ChunkCoord> colliderUnion = Union(colliderRequests);
            HashSet<ChunkCoord> wanted = new HashSet<ChunkCoord>(renderUnion);
            wanted.UnionWith(colliderUnion);

            // Drop presences nobody wants anymore.
            foreach (ChunkCoord coord in presences.Keys.Where(c => !wanted.Contains(c)).ToList())
                DestroyPresence(coord);

            // Create/update for wanted coords.
            foreach (ChunkCoord coord in wanted)
                EnsurePresence(coord, renderUnion.Contains(coord), colliderUnion.Contains(coord));
        }

        void EnsurePresence(ChunkCoord coord, bool wantRender, bool wantCollider)
        {
            PrimalChunk primal = LookupPrimal(coord);
            if (primal == null) return; // not arrived in any mirror yet; try again next Apply.

            // (a) Make sure the chunk has a cached dual at the current version. The dual
            // belongs to the chunk, not the presence — it survives presentation churn and is
            // invalidated only by version changes (own or neighbour, via TerrainModel.Install).
            if (primal.Dual == null || primal.DualBuiltFromVersion != primal.Version)
            {
                BuildNeighborBorderMaps(coord, out var neighborFaces, out var minNeighborCoord);
                var (dualPolygons, _) = PolygonGridGenerator.GenerateDual(
                    primal, hexRadius, neighborFaces, minNeighborCoord);
                primal.Dual = dualPolygons.ToList();
                primal.DualBuiltFromVersion = primal.Version;
            }

            // (b) Make sure a Presence exists with a MeshFilter (always needed).
            if (!presences.TryGetValue(coord, out Presence p))
            {
                GameObject go = new GameObject($"Chunk {coord}");
                go.transform.SetParent(parent, false);
                p = new Presence
                {
                    Go = go,
                    Mf = go.AddComponent<MeshFilter>(),
                };
                presences[coord] = p;
            }

            // (c) Rebuild the Unity Mesh only when the underlying dual has actually changed.
            // Identity check on DualRef catches both the version bump (new list assigned) and
            // the cascade invalidation (Dual replaced with a fresh list).
            if (p.MeshBuiltFromVersion != primal.Version || !ReferenceEquals(p.DualRef, primal.Dual))
            {
                if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
                p.Mesh = BuildMesh(primal.Dual);
                p.Mf.sharedMesh = p.Mesh;
                if (p.Mc != null) p.Mc.sharedMesh = p.Mesh; // collider needs to re-cook
                p.MeshBuiltFromVersion = primal.Version;
                p.DualRef = primal.Dual;
            }

            // (d) Toggle the MeshRenderer.
            if (wantRender && p.Mr == null)
            {
                p.Mr = p.Go.AddComponent<MeshRenderer>();
                p.Mr.sharedMaterial = material;
            }
            else if (!wantRender && p.Mr != null)
            {
                UnityEngine.Object.Destroy(p.Mr);
                p.Mr = null;
            }

            // (e) Toggle the MeshCollider.
            if (wantCollider && p.Mc == null)
            {
                p.Mc = p.Go.AddComponent<MeshCollider>();
                p.Mc.sharedMesh = p.Mesh;
            }
            else if (!wantCollider && p.Mc != null)
            {
                UnityEngine.Object.Destroy(p.Mc);
                p.Mc = null;
            }
        }

        void DestroyPresence(ChunkCoord coord)
        {
            if (!presences.TryGetValue(coord, out Presence p)) return;
            if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
            UnityEngine.Object.Destroy(p.Go);
            presences.Remove(coord);
        }

        public void Clear()
        {
            foreach (Presence p in presences.Values)
            {
                if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
                UnityEngine.Object.Destroy(p.Go);
            }
            presences.Clear();
            renderRequests.Clear();
            colliderRequests.Clear();
        }

        // --- helpers ---

        PrimalChunk LookupPrimal(ChunkCoord coord)
        {
            for (int i = 0; i < primalSources.Count; i++)
            {
                PrimalChunk pc = primalSources[i](coord);
                if (pc != null) return pc;
            }
            return null;
        }

        static HashSet<ChunkCoord> Union(Dictionary<object, HashSet<ChunkCoord>> requests)
        {
            HashSet<ChunkCoord> result = new HashSet<ChunkCoord>();
            foreach (HashSet<ChunkCoord> set in requests.Values) result.UnionWith(set);
            return result;
        }

        // Gather, from the chunk's 6 loaded neighbours (any primal source), the faces
        // incident to each shared border vertex (for cell completion) and the lowest
        // neighbour coord per border position (for deterministic ownership — the chunk
        // with the smallest coord among contributors owns the cell, independent of which
        // chunks are presented). Keyed by PolygonGridGenerator.LatticeKey.
        void BuildNeighborBorderMaps(ChunkCoord coord,
            out Dictionary<(int, int), List<Polygon>> neighborFaces,
            out Dictionary<(int, int), ChunkCoord> minNeighborCoord)
        {
            neighborFaces = new Dictionary<(int, int), List<Polygon>>();
            minNeighborCoord = new Dictionary<(int, int), ChunkCoord>();

            foreach (ChunkCoord n in coord.HexesInRange(1))
            {
                if (n == coord) continue;
                PrimalChunk np = LookupPrimal(n);
                if (np == null) continue;

                foreach (Vertex v in np.Verts.ToArray())
                {
                    if (!v.IsEdge) continue;
                    (int, int) key = PolygonGridGenerator.LatticeKey(v.Position, hexRadius);

                    if (!neighborFaces.TryGetValue(key, out List<Polygon> list))
                        neighborFaces[key] = list = new List<Polygon>();
                    list.AddRange(v.Polygons);

                    if (!minNeighborCoord.TryGetValue(key, out ChunkCoord cur) || n.CompareTo(cur) < 0)
                        minNeighborCoord[key] = n;
                }
            }
        }

        static Mesh BuildMesh(List<Polygon> polygons)
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
            return mesh;
        }
    }
}
