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
    // Catalog-driven overlays: each TileKind's visual is owned by a TilePresenter
    // ScriptableObject (RoadPresenter, BuildingPresenter, etc., resolved through the
    // TileCatalog). The surface iterates the catalog per Apply to rebuild stale overlays
    // and tear down overlays whose chunk has dropped to zero cells of that kind. Adding a
    // new tile kind = author a new Presenter + TileDefinition asset; no edits here.
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

            // Sparse per-kind overlay state. Sized lazily — only kinds that have ever
            // had a cell in this chunk get a Handle entry. Owned by ChunkSurface;
            // presenters mutate handle fields but don't keep their own references.
            public readonly Dictionary<TileKind, ChunkOverlayHandle> Overlays =
                new Dictionary<TileKind, ChunkOverlayHandle>();

            // Outline overlay: a child GameObject under `Go` carrying a line-topology
            // Mesh that traces every dual cell's wedge perimeter via V → M segments
            // per dual polygon edge. Not a per-tile overlay — it's a global RTS-mode
            // decorator over the whole terrain, so it stays separate from the catalog
            // path. Lazily instantiated, version-gated, visibility tracks ShowOutlines.
            public GameObject OutlineGo;
            public MeshFilter OutlineMf;
            public MeshRenderer OutlineMr;
            public Mesh OutlineMesh;
            public int OutlineBuiltFromVersion = -1;
        }

        readonly Transform parent;
        readonly Material material;
        readonly Material outlineMaterial;
        readonly TileCatalog catalog;
        readonly float hexRadius;
        readonly int chunkGridSize;
        readonly int presenceCacheSize;

        // Global outline visibility. GameModeController flips this on Tab; the
        // next Apply() pushes it into every per-chunk OutlineMr.enabled. Public
        // so the controller can also call RequestOutlineRefresh() to make the
        // change visible on the same frame.
        public static bool ShowOutlines = false;

        // Reused scratch MaterialPropertyBlock for setting per-chunk shader
        // uniforms (_ChunkQR) on each terrain MeshRenderer.
        readonly MaterialPropertyBlock chunkMpb = new MaterialPropertyBlock();
        static readonly int ChunkQRId = Shader.PropertyToID("_ChunkQR");

        readonly List<Func<ChunkCoord, PrimalChunk>> primalSources = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> renderRequests = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> colliderRequests = new();
        readonly Dictionary<ChunkCoord, Presence> presences = new();

        // Warm cache: inactive GameObjects whose underlying chunk left the wanted set.
        readonly Dictionary<ChunkCoord, Presence> cache = new();
        readonly Queue<ChunkCoord> cacheOrder = new();

        Mesh flatTile;

        public ChunkSurface(Transform parent, Material material, Material outlineMaterial,
                            TileCatalog catalog,
                            float hexRadius, int chunkGridSize,
                            int presenceCacheSize = 64)
        {
            this.parent = parent;
            this.material = material;
            this.outlineMaterial = outlineMaterial;
            this.catalog = catalog;
            this.hexRadius = hexRadius;
            this.chunkGridSize = chunkGridSize;
            this.presenceCacheSize = Mathf.Max(0, presenceCacheSize);
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

        public void SetRenderRequest(object client, IEnumerable<ChunkCoord> coords)
        {
            if (coords == null) renderRequests.Remove(client);
            else renderRequests[client] = new HashSet<ChunkCoord>(coords);
        }

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

            using (TerrainProfiler.Measure(TerrainProfiler.Phase.ApplyTrim))
                TrimCache();

            using (TerrainProfiler.Measure(TerrainProfiler.Phase.ApplyCount))
            {
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
            if (!primal.IsFlat && primal.Dual == null) return;

            if (!presences.TryGetValue(coord, out Presence p))
            {
                if (cache.TryGetValue(coord, out p))
                {
                    cache.Remove(coord);
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

            // ---- Terrain mesh (flat fast path vs. per-chunk dual mesh) ----
            Mesh activeMesh;
            if (primal.IsFlat)
            {
                if (!p.Flat)
                {
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
                    p.Go.transform.localPosition = Vector3.zero;
                    p.Mf.sharedMesh = null;
                    p.Flat = false;
                    p.MeshBuiltFromVersion = -1;
                }

                if (p.MeshBuiltFromVersion != primal.Version)
                {
                    if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
                    using (TerrainProfiler.Measure(TerrainProfiler.Phase.BuildMesh))
                        p.Mesh = BuildMesh(primal);
                    using (TerrainProfiler.Measure(TerrainProfiler.Phase.MeshAssign))
                        p.Mf.sharedMesh = p.Mesh;
                    if (p.Mc != null)
                        using (TerrainProfiler.Measure(TerrainProfiler.Phase.ColliderCook))
                            p.Mc.sharedMesh = p.Mesh;
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
                // Stamp the chunk's axial coord onto the renderer via MPB so the
                // shader can compare against the hovered chunk and only tint cells
                // in the correct chunk. Set once at attach time — coord doesn't
                // change for the renderer's lifetime. Breaks SRP batching for
                // terrain (acceptable: ~50 visible chunks).
                chunkMpb.Clear();
                chunkMpb.SetVector(ChunkQRId, new Vector4(coord.q, coord.r, 0f, 0f));
                p.Mr.SetPropertyBlock(chunkMpb);
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

            // ---- Per-kind overlays (catalog-driven) ----
            //
            // Flat chunks skip — they mount the shared flat tile, can't host paintable
            // land tiles, and don't ship a primal neighbour table to read from. We also
            // tear down every overlay when wantRender drops (collider-only consumer),
            // matching the prior Road/Building behaviour.
            if (catalog != null && !primal.IsFlat)
            {
                if (wantRender)
                {
                    foreach (TileDefinition def in catalog.Definitions)
                    {
                        if (def == null || def.presenter == null) continue;

                        if (!p.Overlays.TryGetValue(def.kind, out ChunkOverlayHandle handle))
                        {
                            handle = new ChunkOverlayHandle();
                            p.Overlays[def.kind] = handle;
                        }

                        PresenterContext ctx = new PresenterContext(
                            coord, primal, def.kind, p.Go.transform,
                            def.material, handle, LookupPrimal,
                            $"{def.displayName}Overlay");
                        bool hasAny = def.presenter.Rebuild(in ctx);

                        // Collider lifecycle: presenters that report RequiresCollider get a
                        // MeshCollider attached on the overlay GO when both wantCollider AND
                        // hasAny hold. Tear down otherwise — the prior Road/Building code did
                        // the same dance.
                        if (def.presenter.RequiresCollider && hasAny && wantCollider
                            && handle.Go != null && handle.Mesh != null)
                        {
                            if (handle.Mc == null)
                            {
                                handle.Mc = handle.Go.AddComponent<MeshCollider>();
                                using (TerrainProfiler.Measure(TerrainProfiler.Phase.ColliderCook))
                                    handle.Mc.sharedMesh = handle.Mesh;
                                TerrainProfiler.IncColliderAssigns();
                            }
                            else if (handle.Mc.sharedMesh != handle.Mesh)
                            {
                                using (TerrainProfiler.Measure(TerrainProfiler.Phase.ColliderCook))
                                    handle.Mc.sharedMesh = handle.Mesh;
                            }
                        }
                        else if (handle.Mc != null)
                        {
                            UnityEngine.Object.Destroy(handle.Mc);
                            handle.Mc = null;
                        }
                    }
                }
                else
                {
                    // Collider-only / off-screen: tear down every overlay's GO so the warm
                    // cache doesn't pin GPU meshes for chunks that aren't drawn anyway.
                    foreach (KeyValuePair<TileKind, ChunkOverlayHandle> kv in p.Overlays)
                    {
                        TileDefinition def = catalog.Get(kv.Key);
                        def?.presenter?.DestroyOverlay(kv.Value);
                    }
                }
            }

            // ---- Outline overlay ----
            if (!primal.IsFlat && wantRender)
            {
                if (p.OutlineGo == null)
                {
                    p.OutlineGo = new GameObject("OutlineOverlay");
                    p.OutlineGo.transform.SetParent(p.Go.transform, false);
                    p.OutlineMf = p.OutlineGo.AddComponent<MeshFilter>();
                    p.OutlineMr = p.OutlineGo.AddComponent<MeshRenderer>();
                    p.OutlineMr.sharedMaterial = outlineMaterial;
                    p.OutlineMr.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
                    p.OutlineMr.receiveShadows = false;
                    p.OutlineMesh = new Mesh
                    {
                        name = $"OutlineMesh {coord}",
                        indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
                    };
                    p.OutlineMf.sharedMesh = p.OutlineMesh;
                }

                if (p.OutlineBuiltFromVersion != primal.Version)
                {
                    BuildOutlineMesh(primal, p.OutlineMesh);
                    p.OutlineBuiltFromVersion = primal.Version;
                }

                p.OutlineMr.enabled = ShowOutlines;
            }
            else if (!wantRender && p.OutlineGo != null)
            {
                DestroyOutlineOverlay(p);
            }
        }

        void ParkPresence(ChunkCoord coord)
        {
            if (!presences.TryGetValue(coord, out Presence p)) return;
            presences.Remove(coord);

            if (presenceCacheSize <= 0)
            {
                DestroyPresenceObjects(p);
                return;
            }

            if (cache.TryGetValue(coord, out Presence existing))
            {
                cache.Remove(coord);
                DestroyPresenceObjects(existing);
            }

            p.Go.SetActive(false);
            cache[coord] = p;
            cacheOrder.Enqueue(coord);
        }

        void TrimCache()
        {
            while (cache.Count > presenceCacheSize && cacheOrder.Count > 0)
            {
                ChunkCoord coord = cacheOrder.Dequeue();
                if (!cache.TryGetValue(coord, out Presence p)) continue;
                cache.Remove(coord);
                DestroyPresenceObjects(p);
            }
        }

        // Tear down the Unity objects backing a presence. Mesh assets survive
        // GameObject destruction, so each Mesh is destroyed explicitly. Catalog
        // overlays are routed through each presenter's DestroyOverlay so a
        // future presenter with non-trivial cleanup gets called correctly.
        void DestroyPresenceObjects(Presence p)
        {
            if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
            if (p.OutlineMesh != null) UnityEngine.Object.Destroy(p.OutlineMesh);
            if (catalog != null)
            {
                foreach (KeyValuePair<TileKind, ChunkOverlayHandle> kv in p.Overlays)
                {
                    TileDefinition def = catalog.Get(kv.Key);
                    if (def?.presenter != null) def.presenter.DestroyOverlay(kv.Value);
                    else
                    {
                        // Catalog changed mid-session; fall back to manual cleanup so the
                        // mesh asset doesn't leak.
                        if (kv.Value.Mesh != null) UnityEngine.Object.Destroy(kv.Value.Mesh);
                        if (kv.Value.Go != null) UnityEngine.Object.Destroy(kv.Value.Go);
                    }
                }
            }
            p.Overlays.Clear();
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
        // once, reused by every flat chunk's GameObject.
        Mesh GetFlatTile()
        {
            if (flatTile != null) return flatTile;
            float R = hexRadius * chunkGridSize;
            float h = Mathf.Sqrt(3f) / 2f;
            Vector3[] verts =
            {
                new Vector3( R,        0f,  0f      ),
                new Vector3( R * 0.5f, 0f,  R * h   ),
                new Vector3(-R * 0.5f, 0f,  R * h   ),
                new Vector3(-R,        0f,  0f      ),
                new Vector3(-R * 0.5f, 0f, -R * h   ),
                new Vector3( R * 0.5f, 0f, -R * h   ),
            };
            int[] tris = { 0, 2, 1,  0, 3, 2,  0, 4, 3,  0, 5, 4 };
            flatTile = new Mesh { name = "ChunkFlatTile", vertices = verts, triangles = tris };
            flatTile.RecalculateNormals();
            flatTile.RecalculateBounds();
            return flatTile;
        }

        // ----------------- Terrain mesh (per-corner wedge decomposition) -----------------
        //
        // Walks each dual polygon and emits the wedge of every non-carved corner. A corner
        // is "carved" if its owning primal cell's kind has TileDefinition.carvesTerrain
        // set — in which case the kind's presenter paints the wedge instead. Carved
        // wedges (Road / Building / Field today) are skipped; non-carved wedges
        // (Default + Market / Lighthouse / Dock) are emitted with a neutral white tint
        // so the structure rendered on top sits over visible ground.
        //
        // Vertices are unshared per triangle so the shader can compute nearest-corner
        // barycentric tinting from the (cornerA, cornerB, cornerC) triple. Today every
        // wedge uses a degenerate triple (all same colour, all same primal index) since
        // we paint the whole wedge in one hue — leaves room for cross-corner gradients
        // if we ever want them.
        Mesh BuildMesh(PrimalChunk primal)
        {
            List<Polygon> polygons = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;

            var vertices = new List<Vector3>();
            var triangles = new List<int>();
            var colors = new List<Color32>();
            var barycentric = new List<Vector3>();
            var cornerA_rgb = new List<Vector3>();
            var cornerB_rgb = new List<Vector3>();
            var cornerC_rgb = new List<Vector3>();
            var cornerPrimal = new List<Vector3>();

            // The terrain wedge tint is currently constant — structures paint themselves
            // and area presenters carve, so there's no need for per-kind base colour.
            // Reserved for future "tint terrain under a structure" extensions.
            Color32 neutral = new Color32(255, 255, 255, 255);
            Vector3 neutralRgb = new Vector3(1f, 1f, 1f);

            int cornerCursor = 0;
            foreach (Polygon poly in polygons)
            {
                Vector3[] verts = poly.GetVerticesPosition();
                if (verts.Length < 3) { cornerCursor += verts.Length; continue; }

                int n = verts.Length;
                int[] cornerPrimalIdx = new int[n];
                bool[] cornerIsCarved = new bool[n];
                for (int i = 0; i < n; i++)
                {
                    int primalIdx = dvpi != null && cornerCursor + i < dvpi.Length
                        ? dvpi[cornerCursor + i] : -1;
                    cornerPrimalIdx[i] = primalIdx;
                    bool carved = false;
                    if (primalIdx >= 0 && tp != null && primalIdx < tp.Length && catalog != null)
                        carved = catalog.IsCarvingKind(tp[primalIdx].Kind);
                    cornerIsCarved[i] = carved;
                }
                cornerCursor += n;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += verts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (cornerIsCarved[i]) continue;

                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 Mprev = 0.5f * (verts[iPrev] + verts[i]);
                    Vector3 Mnext = 0.5f * (verts[i] + verts[iNext]);

                    Vector3 baryA = new Vector3(1, 0, 0);
                    Vector3 baryB = new Vector3(0, 1, 0);
                    Vector3 baryC = new Vector3(0, 0, 1);
                    int idx = cornerPrimalIdx[i];
                    Vector3 vP = new Vector3(idx, idx, idx);

                    int baseIdx = vertices.Count;
                    vertices.Add(verts[i]); colors.Add(neutral); barycentric.Add(baryA);
                    cornerA_rgb.Add(neutralRgb); cornerB_rgb.Add(neutralRgb); cornerC_rgb.Add(neutralRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(Mnext);    colors.Add(neutral); barycentric.Add(baryB);
                    cornerA_rgb.Add(neutralRgb); cornerB_rgb.Add(neutralRgb); cornerC_rgb.Add(neutralRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(V);        colors.Add(neutral); barycentric.Add(baryC);
                    cornerA_rgb.Add(neutralRgb); cornerB_rgb.Add(neutralRgb); cornerC_rgb.Add(neutralRgb);
                    cornerPrimal.Add(vP);
                    triangles.Add(baseIdx);
                    triangles.Add(baseIdx + 1);
                    triangles.Add(baseIdx + 2);

                    int baseIdx2 = vertices.Count;
                    vertices.Add(verts[i]); colors.Add(neutral); barycentric.Add(baryA);
                    cornerA_rgb.Add(neutralRgb); cornerB_rgb.Add(neutralRgb); cornerC_rgb.Add(neutralRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(Mprev);    colors.Add(neutral); barycentric.Add(baryB);
                    cornerA_rgb.Add(neutralRgb); cornerB_rgb.Add(neutralRgb); cornerC_rgb.Add(neutralRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(V);        colors.Add(neutral); barycentric.Add(baryC);
                    cornerA_rgb.Add(neutralRgb); cornerB_rgb.Add(neutralRgb); cornerC_rgb.Add(neutralRgb);
                    cornerPrimal.Add(vP);
                    triangles.Add(baseIdx2);
                    triangles.Add(baseIdx2 + 1);
                    triangles.Add(baseIdx2 + 2);
                }
            }

            Mesh mesh = new Mesh
            {
                indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
                vertices = vertices.ToArray(),
                triangles = triangles.ToArray(),
                colors32 = colors.ToArray(),
            };
            mesh.SetUVs(0, barycentric);
            mesh.SetUVs(1, cornerA_rgb);
            mesh.SetUVs(2, cornerB_rgb);
            mesh.SetUVs(3, cornerC_rgb);
            mesh.SetUVs(4, cornerPrimal);
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();
            return mesh;
        }

        // ----------------- Outline overlay (wedge-perimeter line topology) ----------
        //
        // For every dual polygon, emits one line segment per polygon edge V → M, where
        // V is the polygon centroid (≈ primal vertex this dual cell surrounds) and M
        // is the midpoint of that polygon edge (= midpoint of two adjacent corner
        // primal centroids). Two adjacent dual polygons each contribute one V → M
        // segment at the shared midpoint M, together forming the full wedge perimeter
        // polyline V₁ → M → V₂ between the two primal cells the edge separates.
        //
        // These segments coincide exactly with the rendered wedge boundaries — every
        // tile presenter is built around the same wedge decomposition.
        static void BuildOutlineMesh(PrimalChunk primal, Mesh mesh)
        {
            List<Polygon> dualPolys = primal.Dual;
            if (dualPolys == null) { mesh.Clear(); return; }

            int segCount = 0;
            foreach (Polygon poly in dualPolys)
            {
                int n = poly.GetVerticesPosition().Length;
                if (n < 3) continue;
                segCount += n;
            }
            if (segCount == 0) { mesh.Clear(); return; }

            Vector3[] verts = new Vector3[segCount * 2];
            int[] idx = new int[segCount * 2];
            Vector3 lift = new Vector3(0f, 0.02f, 0f);

            int vCursor = 0;
            int iCursor = 0;
            foreach (Polygon poly in dualPolys)
            {
                Vector3[] pverts = poly.GetVerticesPosition();
                int n = pverts.Length;
                if (n < 3) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;
                Vector3 Vlift = V + lift;

                for (int i = 0; i < n; i++)
                {
                    int next = (i + 1) % n;
                    Vector3 M = 0.5f * (pverts[i] + pverts[next]) + lift;
                    verts[vCursor] = Vlift;
                    verts[vCursor + 1] = M;
                    idx[iCursor] = vCursor;
                    idx[iCursor + 1] = vCursor + 1;
                    vCursor += 2;
                    iCursor += 2;
                }
            }

            mesh.Clear();
            mesh.vertices = verts;
            mesh.SetIndices(idx, MeshTopology.Lines, 0);
            mesh.RecalculateBounds();
        }

        static void DestroyOutlineOverlay(Presence p)
        {
            if (p.OutlineMesh != null) { UnityEngine.Object.Destroy(p.OutlineMesh); p.OutlineMesh = null; }
            if (p.OutlineGo != null) { UnityEngine.Object.Destroy(p.OutlineGo); p.OutlineGo = null; }
            p.OutlineMf = null;
            p.OutlineMr = null;
            p.OutlineBuiltFromVersion = -1;
        }

        // Push the current ShowOutlines flag into every live presence's outline
        // renderer in O(presences). Used by GameModeController on Tab toggle so
        // the swap is visible on the same frame.
        public void RequestOutlineRefresh()
        {
            foreach (Presence p in presences.Values)
                if (p.OutlineMr != null) p.OutlineMr.enabled = ShowOutlines;
        }
    }
}
