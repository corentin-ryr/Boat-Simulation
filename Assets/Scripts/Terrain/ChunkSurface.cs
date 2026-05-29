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
    // ChunkCoord.MinimalHalo.
    //
    // Phase 2: the dual is built by the streamer worker (in parallel, alongside relaxation)
    // and travels on the PrimalChunk through the publish pipeline. The surface no longer
    // computes any polygon math — it just builds Unity Mesh objects from primal.Dual when
    // the chunk's version changes. All cross-chunk lookups, ownership arbitration, and
    // cascade invalidation happen on the worker thread.
    public class ChunkSurface
    {
        class Presence
        {
            public GameObject Go;
            public MeshFilter Mf;
            public MeshRenderer Mr;     // null when not currently rendered
            public MeshCollider Mc;     // null when not currently collided
            public Mesh Mesh;
            public int MeshBuiltFromVersion = -1;
        }

        readonly Transform parent;
        readonly Material material;

        readonly List<Func<ChunkCoord, PrimalChunk>> primalSources = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> renderRequests = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> colliderRequests = new();
        readonly Dictionary<ChunkCoord, Presence> presences = new();

        public ChunkSurface(Transform parent, Material material, float hexRadius)
        {
            this.parent = parent;
            this.material = material;
            // hexRadius is no longer used by the surface (dual is pre-built by the worker)
            // but kept in the signature to avoid churn at construction sites.
            _ = hexRadius;
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
        // safe to call every frame.
        public void Apply()
        {
            HashSet<ChunkCoord> renderUnion = Union(renderRequests);
            HashSet<ChunkCoord> colliderUnion = Union(colliderRequests);
            HashSet<ChunkCoord> wanted = new HashSet<ChunkCoord>(renderUnion);
            wanted.UnionWith(colliderUnion);

            foreach (ChunkCoord coord in presences.Keys.Where(c => !wanted.Contains(c)).ToList())
                DestroyPresence(coord);

            foreach (ChunkCoord coord in wanted)
                EnsurePresence(coord, renderUnion.Contains(coord), colliderUnion.Contains(coord));
        }

        void EnsurePresence(ChunkCoord coord, bool wantRender, bool wantCollider)
        {
            PrimalChunk primal = LookupPrimal(coord);
            // Either the primal hasn't arrived yet, or it arrived but the worker hasn't
            // built its dual yet (rare — happens briefly during initial streaming bursts).
            // Skip; we'll try again next Apply.
            if (primal == null || primal.Dual == null) return;

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

            // Rebuild the Unity Mesh only when the chunk's version moves. Version covers
            // both primal changes (relaxation moved vertices) and dual changes (a neighbour
            // cascade invalidated our cell ownership) — both bump Version on the worker.
            if (p.MeshBuiltFromVersion != primal.Version)
            {
                if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
                p.Mesh = BuildMesh(primal.Dual);
                p.Mf.sharedMesh = p.Mesh;
                if (p.Mc != null) p.Mc.sharedMesh = p.Mesh; // collider needs to re-cook
                p.MeshBuiltFromVersion = primal.Version;
            }

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
