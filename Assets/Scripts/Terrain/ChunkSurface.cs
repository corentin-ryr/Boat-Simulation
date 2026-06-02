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
                    // Drop the GameObject by SeabedY so the shared tile sits at the seabed level
                    // — the dual mesh path bakes Y per-vertex via Elevation.Sample, but the
                    // flat-tile path uses one shared mesh whose verts are pinned at local Y=0,
                    // so the depth has to come from the transform.
                    if (p.Mesh != null) { UnityEngine.Object.Destroy(p.Mesh); p.Mesh = null; }
                    Vector3 flatCenter = coord.WorldCenter(hexRadius, chunkGridSize);
                    flatCenter.y = Elevation.SeabedY;
                    p.Go.transform.localPosition = flatCenter;
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
                        p.Mesh = BuildMesh(primal);
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

        // Static palette mapping a TileKind to its vertex-colour multiplier. Neutral white
        // is the no-op multiply, so untouched cells render exactly as before. Painted kinds
        // shift the base hue without dimming it too aggressively (alpha=255 for all). The
        // sentinel entry at index `kindCount` covers neighbour-chunk corners (primal index
        // -1) — same neutral as Default so unowned border vertices blend invisibly.
        static readonly Color32[] TileKindPalette =
        {
            new Color32(255, 255, 255, 255), // Default — neutral, base colour unchanged
            new Color32(160, 150, 135, 255), // Road    — warm grey
            new Color32(210, 120, 110, 255), // Building — soft red
            new Color32(255, 255, 255, 255), // [-1 sentinel] — no owning primal, neutral
        };
        static readonly int PaletteNeutralIndex = TileKindPalette.Length - 1;

        static Color32 ColorFor(int primalIdx, TileProperty[] tp)
        {
            if (primalIdx < 0 || tp == null || primalIdx >= tp.Length)
                return TileKindPalette[PaletteNeutralIndex];
            int kind = (int)tp[primalIdx].Kind;
            if ((uint)kind >= (uint)PaletteNeutralIndex) return TileKindPalette[PaletteNeutralIndex];
            return TileKindPalette[kind];
        }

        // RGB components of a palette entry, packed for upload as a UV channel. Alpha is
        // dropped since the shader only multiplies RGB; this also keeps the channel as a
        // Vector3 instead of Vector4.
        static Vector3 ColorRgb(Color32 c) => new Vector3(c.r / 255f, c.g / 255f, c.b / 255f);

        static Mesh BuildMesh(PrimalChunk primal)
        {
            List<Polygon> polygons = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;          // corner → primal cell index (or -1)
            TileProperty[] tp = primal.TileProperties;          // null = all Default → all neutral

            // We unshare vertices per TRIANGLE so the fragment shader can compute nearest
            // neighbour across the triangle's three primal corners. Each emitted vertex
            // carries the same (cA, cB, cC) triple of corner colours plus a barycentric
            // identifier (which of the three I am). After linear interpolation in the
            // pipeline, the fragment sees a constant (cA, cB, cC) and a per-pixel barycentric
            // weighting; it picks whichever corner has the largest weight. Result: each
            // primal-Voronoi region of the triangle takes its owning cell's palette colour
            // with a sharp edge along the bisectors.
            //
            // Storage cost: a polygon with N corners has (N-2) triangles, each with 3 unique
            // vertices → 3·(N-2) verts where the shared-fan version had N. For typical
            // dual quads (N=4) this is 1.5× the vertex count; for hexagonal-ish duals (N=6)
            // it's 2×. Triangle and index counts are unchanged.
            var vertices    = new List<Vector3>();
            var triangles   = new List<int>();
            var colors      = new List<Color32>();   // per-vertex own colour (back-compat / fallback)
            var barycentric = new List<Vector3>();   // (1,0,0), (0,1,0), (0,0,1) per triangle vertex
            var cornerA_rgb = new List<Vector3>();   // all three identical per triangle, so each
            var cornerB_rgb = new List<Vector3>();   // vertex in the triangle outputs the same triple
            var cornerC_rgb = new List<Vector3>();

            int cornerCursor = 0;   // index into dvpi; advances by each polygon's corner count
            foreach (Polygon p in polygons)
            {
                Vector3[] verts = p.GetVerticesPosition();
                if (verts.Length < 3)
                {
                    // GenerateDual already skipped these, but if a legacy/restored dual still
                    // contains one we must keep dvpi's cursor in sync with the corner walk.
                    cornerCursor += verts.Length;
                    continue;
                }

                // Resolve every corner's palette colour up front so each fan triangle below
                // can pick three of them without redoing the dvpi/tp lookups.
                int n = verts.Length;
                Color32[] cornerColors = new Color32[n];
                for (int i = 0; i < n; i++)
                {
                    int primalIdx = dvpi != null && cornerCursor + i < dvpi.Length
                        ? dvpi[cornerCursor + i] : -1;
                    cornerColors[i] = ColorFor(primalIdx, tp);
                }
                cornerCursor += n;

                // Triangle fan from corner 0: indices (0, i+1, i) for i = 1..n-2. We emit
                // three NEW vertices per triangle; never sharing means each triangle's
                // barycentric identification is self-contained.
                for (int i = 1; i < n - 1; i++)
                {
                    int idxA = 0;
                    int idxB = i + 1;
                    int idxC = i;

                    Color32 cA = cornerColors[idxA];
                    Color32 cB = cornerColors[idxB];
                    Color32 cC = cornerColors[idxC];
                    Vector3 vA = ColorRgb(cA);
                    Vector3 vB = ColorRgb(cB);
                    Vector3 vC = ColorRgb(cC);

                    int baseIdx = vertices.Count;

                    // Vertex A — owns barycentric (1,0,0)
                    vertices.Add(verts[idxA]);
                    colors.Add(cA);
                    barycentric.Add(new Vector3(1, 0, 0));
                    cornerA_rgb.Add(vA); cornerB_rgb.Add(vB); cornerC_rgb.Add(vC);

                    // Vertex B — (0,1,0)
                    vertices.Add(verts[idxB]);
                    colors.Add(cB);
                    barycentric.Add(new Vector3(0, 1, 0));
                    cornerA_rgb.Add(vA); cornerB_rgb.Add(vB); cornerC_rgb.Add(vC);

                    // Vertex C — (0,0,1)
                    vertices.Add(verts[idxC]);
                    colors.Add(cC);
                    barycentric.Add(new Vector3(0, 0, 1));
                    cornerA_rgb.Add(vA); cornerB_rgb.Add(vB); cornerC_rgb.Add(vC);

                    // Same winding the legacy fan used: (0, i+1, i).
                    triangles.Add(baseIdx);
                    triangles.Add(baseIdx + 1);
                    triangles.Add(baseIdx + 2);
                }
            }

            Mesh mesh = new Mesh
            {
                indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
                vertices = vertices.ToArray(),
                triangles = triangles.ToArray(),
                colors32 = colors.ToArray(),
            };
            mesh.SetUVs(0, barycentric);     // TEXCOORD0 → fragment barycentric
            mesh.SetUVs(1, cornerA_rgb);     // TEXCOORD1 → triangle corner A colour
            mesh.SetUVs(2, cornerB_rgb);     // TEXCOORD2 → triangle corner B colour
            mesh.SetUVs(3, cornerC_rgb);     // TEXCOORD3 → triangle corner C colour
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();
            return mesh;
        }
    }
}
