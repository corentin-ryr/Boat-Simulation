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
    // Flat-ocean fast path: when primal.IsFlat is true, the surface mounts a single shared
    // hex-tile mesh (size = chunk circumradius) at the chunk's world center, used for both
    // renderer and collider, instead of building per-cell geometry from the dual. The same
    // Mesh object is reused across every flat chunk — one allocation, one MeshCollider cook,
    // shared by potentially thousands of distant-ocean chunks.
    //
    // Warm-presence cache: when a chunk leaves the wanted set, its Presence (GameObject +
    // Mesh + MeshCollider) is SetActive(false) and parked in a FIFO cache instead of being
    // destroyed. On reentry the GameObject is SetActive(true) and — if primal.Version still
    // matches — the existing Mesh and cooked MeshCollider are reused as-is, skipping the
    // ~10ms-per-chunk BuildMesh + ColliderCook hitch that would otherwise hit the main
    // thread on every camera-pan revisit. Overflow evicts FIFO; size bound is
    // `presenceCacheSize` (constructor argument).
    //
    // Phase 2: for non-flat chunks the dual is built by the streamer worker (in parallel,
    // alongside relaxation) and travels on the PrimalChunk through the publish pipeline. The
    // surface no longer computes any polygon math — it just builds Unity Mesh objects from
    // primal.Dual when the chunk's version changes. All cross-chunk lookups, ownership
    // arbitration, and cascade invalidation happen on the worker thread.
    public class ChunkSurface
    {
        class Presence
        {
            public GameObject Go;
            public MeshFilter Mf;
            public MeshRenderer Mr;     // null when not currently rendered
            public MeshCollider Mc;     // null when not currently collided
            public Mesh Mesh;           // owned per-presence mesh; null when using the shared flat tile
            public int MeshBuiltFromVersion = -1;
            public bool Flat;           // mode flag: true while sharedMesh == flatTile and GO is positioned
        }

        readonly Transform parent;
        readonly Material material;
        readonly float hexRadius;
        readonly int chunkGridSize;
        readonly int presenceCacheSize;

        readonly List<Func<ChunkCoord, PrimalChunk>> primalSources = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> renderRequests = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> colliderRequests = new();
        readonly Dictionary<ChunkCoord, Presence> presences = new();

        // Warm cache: inactive GameObjects whose underlying chunk left the wanted set.
        // FIFO eviction (`cacheOrder` may contain stale entries — promoted-out coords — that
        // are skipped during TrimCache).
        readonly Dictionary<ChunkCoord, Presence> cache = new();
        readonly Queue<ChunkCoord> cacheOrder = new();

        Mesh flatTile; // lazy — one shared Mesh covering a single chunk's hex footprint at Y=0

        public ChunkSurface(Transform parent, Material material, float hexRadius, int chunkGridSize,
                            int presenceCacheSize = 64)
        {
            this.parent = parent;
            this.material = material;
            this.hexRadius = hexRadius;
            this.chunkGridSize = chunkGridSize;
            this.presenceCacheSize = Mathf.Max(0, presenceCacheSize);
        }

        // For debug gizmos: the cached duals of every presented chunk's primal. Flat chunks
        // expose their empty sentinel Dual, which yields nothing visible — correct.
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
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.SurfaceApply);
            TerrainProfiler.IncSurfaceApplies();

            HashSet<ChunkCoord> renderUnion = Union(renderRequests);
            HashSet<ChunkCoord> colliderUnion = Union(colliderRequests);
            HashSet<ChunkCoord> wanted = new HashSet<ChunkCoord>(renderUnion);
            wanted.UnionWith(colliderUnion);

            // Move now-unwanted presences into the warm cache (or destroy if the cache is
            // disabled). They stay in the scene as inactive GameObjects, ready for a free
            // revival on the next camera pan that brings them back into view.
            using (TerrainProfiler.Measure(TerrainProfiler.Phase.ApplyPark))
            {
                foreach (ChunkCoord coord in presences.Keys.Where(c => !wanted.Contains(c)).ToList())
                    ParkPresence(coord);
            }

            using (TerrainProfiler.Measure(TerrainProfiler.Phase.ApplyEnsure))
            {
                foreach (ChunkCoord coord in wanted)
                    EnsurePresence(coord, renderUnion.Contains(coord), colliderUnion.Contains(coord));
            }

            // FIFO trim of the cache to the configured bound. Cheap when under threshold;
            // skips stale entries that were already promoted out.
            using (TerrainProfiler.Measure(TerrainProfiler.Phase.ApplyTrim))
                TrimCache();

            using (TerrainProfiler.Measure(TerrainProfiler.Phase.ApplyCount))
            {
                // Live counts for the stats overlay. GetIndexCount is allocation-free (unlike
                // `mesh.triangles.Length`, which materializes a fresh int[] every call).
                TerrainProfiler.Presences = presences.Count;
                TerrainProfiler.PresencesCached = cache.Count;
                int liveTris = 0;
                foreach (Presence p in presences.Values)
                {
                    Mesh m = p.Mf != null ? p.Mf.sharedMesh : null;
                    if (m != null && m.subMeshCount > 0) liveTris += (int)(m.GetIndexCount(0) / 3);
                }
                TerrainProfiler.LiveTriangles = liveTris;
            }
        }

        void EnsurePresence(ChunkCoord coord, bool wantRender, bool wantCollider)
        {
            PrimalChunk primal = LookupPrimal(coord);
            if (primal == null) return;
            // Non-flat chunks need their dual built before we can mesh them; flat chunks use
            // the shared tile and don't depend on Dual at all.
            if (!primal.IsFlat && primal.Dual == null) return;

            if (!presences.TryGetValue(coord, out Presence p))
            {
                // Cache check before fresh allocation: a chunk that was parked while still
                // sharing the current primal Version revives with no mesh work at all — just
                // a SetActive flip and (maybe) toggling renderer/collider components.
                if (cache.TryGetValue(coord, out p))
                {
                    cache.Remove(coord);
                    // stale entry in cacheOrder is harmless; TrimCache will skip it.
                    p.Go.SetActive(true);
                }
                else
                {
                    GameObject go = new GameObject($"Chunk {coord}");
                    go.transform.SetParent(parent, false);
                    p = new Presence
                    {
                        Go = go,
                        Mf = go.AddComponent<MeshFilter>(),
                    };
                }
                presences[coord] = p;
            }

            // Pick the mesh source for this presence. Flat → shared tile (one allocation for
            // the entire world). Non-flat → per-chunk mesh built from the dual, rebuilt when
            // the chunk's Version moves.
            Mesh activeMesh;
            if (primal.IsFlat)
            {
                if (!p.Flat)
                {
                    // Switch to flat mode: release any per-chunk mesh, position the GameObject at
                    // the chunk's world center (the shared tile is in local space around origin).
                    if (p.Mesh != null) { UnityEngine.Object.Destroy(p.Mesh); p.Mesh = null; }
                    p.Go.transform.localPosition = coord.WorldCenter(hexRadius, chunkGridSize);
                    p.Flat = true;
                    p.MeshBuiltFromVersion = primal.Version;
                    p.Mf.sharedMesh = GetFlatTile();
                    if (p.Mc != null) p.Mc.sharedMesh = p.Mf.sharedMesh;
                }
                activeMesh = p.Mf.sharedMesh;
            }
            else
            {
                if (p.Flat)
                {
                    // Switch back to per-chunk mode: reset transform; the dual mesh is in world
                    // coords. Force a rebuild on the next version check.
                    p.Go.transform.localPosition = Vector3.zero;
                    p.Mf.sharedMesh = null;
                    p.Flat = false;
                    p.MeshBuiltFromVersion = -1;
                }

                // Rebuild the Unity Mesh only when the chunk's version moves. Version covers
                // both primal changes (relaxation moved vertices) and dual changes (a neighbour
                // cascade invalidated our cell ownership) — both bump Version on the worker.
                // This is also the fast path for warm-cache revivals: if the cached presence
                // still matches the primal's Version, we skip BuildMesh + ColliderCook entirely.
                if (p.MeshBuiltFromVersion != primal.Version)
                {
                    if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
                    using (TerrainProfiler.Measure(TerrainProfiler.Phase.BuildMesh))
                        p.Mesh = BuildMesh(primal.Dual);
                    using (TerrainProfiler.Measure(TerrainProfiler.Phase.MeshAssign))
                        p.Mf.sharedMesh = p.Mesh;
                    if (p.Mc != null)
                        using (TerrainProfiler.Measure(TerrainProfiler.Phase.ColliderCook))
                            p.Mc.sharedMesh = p.Mesh; // collider needs to re-cook
                    p.MeshBuiltFromVersion = primal.Version;
                    TerrainProfiler.IncMeshRebuilds();
                    if (p.Mesh != null) TerrainProfiler.AddTrianglesBuilt((int)(p.Mesh.GetIndexCount(0) / 3));
                }
                activeMesh = p.Mf.sharedMesh;
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
                using (TerrainProfiler.Measure(TerrainProfiler.Phase.ColliderCook))
                    p.Mc.sharedMesh = activeMesh;
                TerrainProfiler.IncColliderAssigns();
            }
            else if (!wantCollider && p.Mc != null)
            {
                UnityEngine.Object.Destroy(p.Mc);
                p.Mc = null;
            }
        }

        // Park a presence in the warm cache (SetActive(false), keep all components and mesh
        // intact). If the cache is disabled (size 0), destroy immediately — same as before.
        void ParkPresence(ChunkCoord coord)
        {
            if (!presences.TryGetValue(coord, out Presence p)) return;
            presences.Remove(coord);

            if (presenceCacheSize <= 0)
            {
                DestroyPresenceObjects(p);
                return;
            }

            // If a stale cache entry exists for this coord (shouldn't normally happen — coord
            // was in presences — but defensive), free it first.
            if (cache.TryGetValue(coord, out Presence existing))
            {
                cache.Remove(coord);
                DestroyPresenceObjects(existing);
            }

            p.Go.SetActive(false);
            cache[coord] = p;
            cacheOrder.Enqueue(coord);
        }

        // FIFO eviction down to presenceCacheSize. Stale queue entries (promoted-out coords)
        // are skipped, matching the streamer's TrimCache pattern.
        void TrimCache()
        {
            while (cache.Count > presenceCacheSize && cacheOrder.Count > 0)
            {
                ChunkCoord coord = cacheOrder.Dequeue();
                if (!cache.TryGetValue(coord, out Presence p)) continue; // promoted out
                cache.Remove(coord);
                DestroyPresenceObjects(p);
            }
        }

        // Tear down the Unity objects backing a presence. The Mesh (when per-chunk, not the
        // shared flat tile) is a separate Unity asset that survives GameObject destruction,
        // so we Destroy it explicitly.
        static void DestroyPresenceObjects(Presence p)
        {
            if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
            if (p.Go != null) UnityEngine.Object.Destroy(p.Go);
        }

        public void Clear()
        {
            foreach (Presence p in presences.Values) DestroyPresenceObjects(p);
            foreach (Presence p in cache.Values) DestroyPresenceObjects(p);
            presences.Clear();
            cache.Clear();
            cacheOrder.Clear();
            renderRequests.Clear();
            colliderRequests.Clear();
            if (flatTile != null) { UnityEngine.Object.Destroy(flatTile); flatTile = null; }
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

        // Single shared flat-ocean tile: a regular flat-top hexagon at the origin in local
        // space, circumradius = chunkGridSize · hexRadius (matches a chunk's footprint). Built
        // once, reused by every flat chunk's GameObject — one Mesh allocation and one
        // MeshCollider cook for potentially thousands of distant-ocean chunks.
        Mesh GetFlatTile()
        {
            if (flatTile != null) return flatTile;

            float R = hexRadius * chunkGridSize;
            float h = Mathf.Sqrt(3f) / 2f;
            Vector3[] verts =
            {
                new Vector3( R,       0f,  0f      ),
                new Vector3( R * 0.5f, 0f,  R * h  ),
                new Vector3(-R * 0.5f, 0f,  R * h  ),
                new Vector3(-R,       0f,  0f      ),
                new Vector3(-R * 0.5f, 0f, -R * h  ),
                new Vector3( R * 0.5f, 0f, -R * h  ),
            };
            // Triangle fan from v0; winding matches BuildMesh's convention so the front face
            // points up (same as the dual mesh, so the shared material renders identically).
            int[] tris =
            {
                0, 2, 1,
                0, 3, 2,
                0, 4, 3,
                0, 5, 4,
            };

            flatTile = new Mesh { name = "ChunkFlatTile", vertices = verts, triangles = tris };
            flatTile.RecalculateNormals();
            flatTile.RecalculateBounds();
            return flatTile;
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
