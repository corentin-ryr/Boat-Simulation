using System;
using System.Collections.Generic;
using System.Linq;
using TerrainGrid.MeshBuilders;
using UnityEngine;

namespace TerrainGrid
{
    // Owns the main-thread "presence" of terrain chunks in the scene: one GameObject per
    // chunk any client wants, with a MeshRenderer when someone wants it drawn and a
    // MeshCollider when someone wants it physical.
    //
    // Per-kind overlays (Road / Building / Field / Market / Lighthouse / Dock) are
    // dispatched explicitly in EnsurePresence — one call per kind through the shared
    // RebuildOverlay helper. Adding a new TileKind takes three coordinated edits:
    //   1. Add the enum entry in TileProperty.cs.
    //   2. Add a Build*Mesh.cs file under MeshBuilders/.
    //   3. Add one Material field on ChunkManager (passed through TileMaterialSet) +
    //      one RebuildOverlay call here.
    // No ScriptableObjects, no catalog assets — just static dispatch.
    //
    // Flat-ocean fast path: when primal.IsFlat is true, the surface mounts a single shared
    // hex-tile mesh (size = chunk circumradius) at the chunk's world center, used for both
    // renderer and collider, instead of building per-cell geometry from the dual.
    //
    // Warm-presence cache: when a chunk leaves the wanted set, its Presence (GameObject +
    // Mesh + MeshCollider + per-kind handles) is SetActive(false) and parked in a FIFO
    // cache instead of being destroyed. On reentry, if primal.Version still matches, the
    // existing meshes and cooked colliders are reused as-is.
    public class ChunkSurface
    {
        // One bundle so ChunkSurface's constructor + the dispatch helper don't grow
        // a six-Material argument list every time a kind is added.
        public struct TileMaterialSet
        {
            public Material Road;
            public Material Building;
            public Material Field;
            public Material Market;
            public Material Lighthouse;
            public Material Dock;
        }

        class Presence
        {
            public GameObject Go;
            public MeshFilter Mf;
            public MeshRenderer Mr;
            public MeshCollider Mc;
            public Mesh Mesh;
            public int MeshBuiltFromVersion = -1;
            public bool Flat;

            // Per-kind overlay state, keyed by TileKind. Sparse: only kinds that
            // have ever had a cell in this chunk get a handle entry. Same key
            // scheme survives across the catalog rollback so warm-cached presences
            // are reusable.
            public readonly Dictionary<TileKind, OverlayHandle> Overlays =
                new Dictionary<TileKind, OverlayHandle>();

            // Outline overlay: a child GameObject under Go carrying a
            // line-topology Mesh tracing every dual cell's wedge perimeter via
            // V → M segments. Global RTS-mode decorator, not a per-tile overlay
            // — stays inline here. Lazily instantiated, version-gated rebuild;
            // visibility tracks ShowOutlines.
            public GameObject OutlineGo;
            public MeshFilter OutlineMf;
            public MeshRenderer OutlineMr;
            public Mesh OutlineMesh;
            public int OutlineBuiltFromVersion = -1;
        }

        readonly Transform parent;
        readonly Material material;
        readonly Material outlineMaterial;
        readonly TileMaterialSet tileMaterials;
        readonly float hexRadius;
        readonly int chunkGridSize;
        readonly int presenceCacheSize;

        // Global outline visibility. GameModeController flips this on Tab; the
        // next Apply() pushes it into every per-chunk OutlineMr.enabled. Public
        // so the controller can also call RequestOutlineRefresh() to make the
        // change visible on the same frame.
        public static bool ShowOutlines = false;

        readonly MaterialPropertyBlock chunkMpb = new MaterialPropertyBlock();
        static readonly int ChunkQRId = Shader.PropertyToID("_ChunkQR");

        readonly List<Func<ChunkCoord, PrimalChunk>> primalSources = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> renderRequests = new();
        readonly Dictionary<object, HashSet<ChunkCoord>> colliderRequests = new();
        readonly Dictionary<ChunkCoord, Presence> presences = new();

        readonly Dictionary<ChunkCoord, Presence> cache = new();
        readonly Queue<ChunkCoord> cacheOrder = new();

        Mesh flatTile;

        public ChunkSurface(Transform parent, Material material, Material outlineMaterial,
                            TileMaterialSet tileMaterials,
                            float hexRadius, int chunkGridSize,
                            int presenceCacheSize = 64)
        {
            this.parent = parent;
            this.material = material;
            this.outlineMaterial = outlineMaterial;
            this.tileMaterials = tileMaterials;
            this.hexRadius = hexRadius;
            this.chunkGridSize = chunkGridSize;
            this.presenceCacheSize = Mathf.Max(0, presenceCacheSize);
        }

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

            // ---- Terrain base mesh (flat fast path or per-chunk dual mesh) ----
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

            // ---- Per-kind overlays (explicit dispatch, one call per kind) ----
            //
            // Flat chunks skip — they mount the shared flat tile, can't host
            // paintable land tiles, and don't ship a primal neighbour table.
            // wantRender=false (collider-only consumer) destroys every overlay
            // GO so the warm cache doesn't pin GPU meshes for off-screen chunks.
            if (!primal.IsFlat)
            {
                if (wantRender)
                {
                    RebuildOverlay(p, TileKind.Road,       primal, tileMaterials.Road,       true,  wantCollider, RoadMesh.Build,       "RoadOverlay");
                    RebuildOverlay(p, TileKind.Building,   primal, tileMaterials.Building,   true,  wantCollider, BuildingMesh.Build,   "BuildingOverlay");
                    RebuildOverlay(p, TileKind.Field,      primal, tileMaterials.Field,      false, wantCollider, FieldMesh.Build,      "FieldOverlay");
                    RebuildOverlay(p, TileKind.Market,     primal, tileMaterials.Market,     false, wantCollider, MarketMesh.Build,     "MarketOverlay");
                    RebuildOverlay(p, TileKind.Lighthouse, primal, tileMaterials.Lighthouse, true,  wantCollider, LighthouseMesh.Build, "LighthouseOverlay");
                    RebuildOverlay(p, TileKind.Dock,       primal, tileMaterials.Dock,       true,  wantCollider, DockMesh.Build,       "DockOverlay");
                }
                else
                {
                    foreach (KeyValuePair<TileKind, OverlayHandle> kv in p.Overlays)
                        DestroyOverlay(kv.Value);
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

        // Per-kind dispatch helper. Owns the entire lifecycle for one (chunk,
        // TileKind) overlay: lazy-spawn the child GO + MeshFilter +
        // MeshRenderer, version-gate the rebuild, swap the Mesh asset safely,
        // attach/detach the MeshCollider when `requiresCollider && wantCollider`,
        // and destroy the GO when the chunk has zero cells of this kind.
        //
        // The mesh builder is a delegate so this helper doesn't know anything
        // per-kind. Adding a kind = pass a new Build delegate.
        void RebuildOverlay(Presence p, TileKind kind, PrimalChunk primal,
                            Material material,
                            bool requiresCollider, bool wantCollider,
                            Func<PrimalChunk, Func<ChunkCoord, PrimalChunk>, Mesh> build,
                            string label)
        {
            if (!p.Overlays.TryGetValue(kind, out OverlayHandle handle))
            {
                handle = new OverlayHandle();
                p.Overlays[kind] = handle;
            }

            // Skip rebuild when same Version and we already have a mesh — the
            // worker bumps Version on any tile mutation (cross-chunk cascade
            // included), so this catches the dominant "nothing changed" case.
            bool versionMatches = handle.BuiltFromVersion == primal.Version && handle.Mesh != null;

            if (!versionMatches)
            {
                Mesh mesh = build(primal, LookupPrimal);
                bool hasGeometry = mesh != null && mesh.vertexCount > 0;

                if (!hasGeometry)
                {
                    // No cells of this kind in the chunk — tear the GO down
                    // entirely. Frees GPU resources and keeps the scene tree
                    // clean. Free the throwaway empty mesh too.
                    if (mesh != null) UnityEngine.Object.Destroy(mesh);
                    DestroyOverlay(handle);
                    handle.BuiltFromVersion = primal.Version;
                    return;
                }

                if (handle.Go == null)
                {
                    handle.Go = new GameObject(label);
                    handle.Go.transform.SetParent(p.Go.transform, false);
                    handle.Mf = handle.Go.AddComponent<MeshFilter>();
                    handle.Mr = handle.Go.AddComponent<MeshRenderer>();
                    handle.Mr.sharedMaterial = material;
                }

                if (handle.Mesh != null) UnityEngine.Object.Destroy(handle.Mesh);
                handle.Mesh = mesh;
                handle.Mesh.name = $"{label} (v{primal.Version})";
                handle.Mf.sharedMesh = handle.Mesh;
                handle.BuiltFromVersion = primal.Version;
            }

            // Collider lifecycle. Attach when the presenter requires it AND
            // the chunk wants colliders AND the overlay actually has geometry.
            // Detach when any of those drops. Re-cook when the Mesh asset
            // changed (different reference).
            if (requiresCollider && wantCollider && handle.Go != null && handle.Mesh != null)
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

        static void DestroyOverlay(OverlayHandle handle)
        {
            if (handle == null) return;
            if (handle.Mc != null) { UnityEngine.Object.Destroy(handle.Mc); handle.Mc = null; }
            if (handle.Mesh != null) { UnityEngine.Object.Destroy(handle.Mesh); handle.Mesh = null; }
            if (handle.Go != null) { UnityEngine.Object.Destroy(handle.Go); handle.Go = null; }
            handle.Mf = null;
            handle.Mr = null;
            handle.BuiltFromVersion = -1;
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

        static void DestroyPresenceObjects(Presence p)
        {
            if (p.Mesh != null) UnityEngine.Object.Destroy(p.Mesh);
            if (p.OutlineMesh != null) UnityEngine.Object.Destroy(p.OutlineMesh);
            foreach (KeyValuePair<TileKind, OverlayHandle> kv in p.Overlays)
                DestroyOverlay(kv.Value);
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

        // Carve list: cells of these kinds have their terrain wedge SKIPPED in
        // BuildMesh, because the corresponding overlay paints the cell
        // pixel-for-pixel. Default + Market/Lighthouse/Dock do NOT carve —
        // those structures sit on top of natural ground.
        static bool IsCarvingKind(TileKind k)
            => k == TileKind.Road || k == TileKind.Building || k == TileKind.Field;

        // Terrain mesh build. Walks every dual polygon and emits the wedge of
        // each non-carved corner. Each wedge is a 2-triangle fan written with
        // unshared vertices (3 vertices per triangle) so the shader's
        // nearest-corner barycentric tinting can paint sharp per-cell
        // Voronoi regions. cornerA/B/C_rgb (TEXCOORD1/2/3) are read by
        // TerrainElevation.shader as a multiplicative tint on top of the
        // height-banded biome colour.
        //
        // Today every wedge writes the SAME triple into all three corner
        // slots (degenerate NN) — gives one hue per wedge. Future cross-corner
        // gradients would write three different colours per triangle.
        static Mesh BuildMesh(PrimalChunk primal)
        {
            List<Polygon> polygons = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;

            var vertices    = new List<Vector3>();
            var triangles   = new List<int>();
            var colors      = new List<Color32>();
            var barycentric = new List<Vector3>();
            var cornerA_rgb = new List<Vector3>();
            var cornerB_rgb = new List<Vector3>();
            var cornerC_rgb = new List<Vector3>();
            var cornerPrimal = new List<Vector3>();

            int cornerCursor = 0;
            foreach (Polygon poly in polygons)
            {
                Vector3[] verts = poly.GetVerticesPosition();
                if (verts.Length < 3)
                {
                    cornerCursor += verts.Length;
                    continue;
                }

                int n = verts.Length;
                Color32[] cornerColors = new Color32[n];
                int[] cornerPrimalIdx = new int[n];
                bool[] cornerIsCarved = new bool[n];
                for (int i = 0; i < n; i++)
                {
                    int primalIdx = dvpi != null && cornerCursor + i < dvpi.Length
                        ? dvpi[cornerCursor + i] : -1;
                    cornerColors[i] = TileKindPalette.ColorFor(primalIdx, tp);
                    cornerPrimalIdx[i] = primalIdx;
                    bool carved = primalIdx >= 0 && tp != null && primalIdx < tp.Length
                                  && IsCarvingKind(tp[primalIdx].Kind);
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

                    Color32 c = cornerColors[i];
                    Vector3 cRgb = TileKindPalette.ColorRgb(c);
                    Vector3 baryA = new Vector3(1, 0, 0);
                    Vector3 baryB = new Vector3(0, 1, 0);
                    Vector3 baryC = new Vector3(0, 0, 1);
                    int idx = cornerPrimalIdx[i];
                    Vector3 vP = new Vector3(idx, idx, idx);

                    // Triangle 1 — (c_i, V, M_next). Winding chosen to match
                    // RoadMesh / BuildingMesh / FieldMesh (and all the new
                    // wedge structure meshes) so terrain wedges face up
                    // under Unity's left-handed front-face convention.
                    // Triangle 2 below does (c_i, M_prev, V) for the same
                    // reason. Both pair up to cover the full corner-wedge
                    // quad (c_i, M_next, V, M_prev) viewed from above.
                    int baseIdx = vertices.Count;
                    vertices.Add(verts[i]); colors.Add(c); barycentric.Add(baryA);
                    cornerA_rgb.Add(cRgb); cornerB_rgb.Add(cRgb); cornerC_rgb.Add(cRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(Mnext);    colors.Add(c); barycentric.Add(baryB);
                    cornerA_rgb.Add(cRgb); cornerB_rgb.Add(cRgb); cornerC_rgb.Add(cRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(V);        colors.Add(c); barycentric.Add(baryC);
                    cornerA_rgb.Add(cRgb); cornerB_rgb.Add(cRgb); cornerC_rgb.Add(cRgb);
                    cornerPrimal.Add(vP);
                    triangles.Add(baseIdx);
                    triangles.Add(baseIdx + 2);    // V before Mnext flips winding
                    triangles.Add(baseIdx + 1);    // Mnext after V — now (c_i, V, Mnext)

                    int baseIdx2 = vertices.Count;
                    vertices.Add(verts[i]); colors.Add(c); barycentric.Add(baryA);
                    cornerA_rgb.Add(cRgb); cornerB_rgb.Add(cRgb); cornerC_rgb.Add(cRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(Mprev);    colors.Add(c); barycentric.Add(baryB);
                    cornerA_rgb.Add(cRgb); cornerB_rgb.Add(cRgb); cornerC_rgb.Add(cRgb);
                    cornerPrimal.Add(vP);
                    vertices.Add(V);        colors.Add(c); barycentric.Add(baryC);
                    cornerA_rgb.Add(cRgb); cornerB_rgb.Add(cRgb); cornerC_rgb.Add(cRgb);
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

        // Wedge-perimeter outline mesh. One V → M segment per dual polygon
        // edge, where V = polygon centroid and M = midpoint of that polygon
        // edge. Two adjacent polygons each contribute one V → M segment at
        // the shared midpoint, together forming the full wedge perimeter
        // polyline. Coincides exactly with rendered cell boundaries.
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

        public void RequestOutlineRefresh()
        {
            foreach (Presence p in presences.Values)
                if (p.OutlineMr != null) p.OutlineMr.enabled = ShowOutlines;
        }
    }
}
