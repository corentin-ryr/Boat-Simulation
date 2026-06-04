using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    // Streaming orchestrator (a MonoBehaviour). Watches the camera and acts as one consumer of
    // the TerrainStreamer: each time the camera crosses a chunk boundary it declares the render
    // set (renderRadius around the camera) via SetDesired. The streamer (a background thread)
    // loads that set plus a 1-ring halo, relaxes the seams, and publishes ready PrimalChunk
    // copies; the main thread installs them into a render-side *mirror* model and meshes the
    // render set through ChunkRenderer.
    //
    // The streamer is multi-consumer: other systems (e.g. a future simulation MonoBehaviour) can
    // register their own consumers and request chunks independently — a chunk stays loaded while
    // any consumer wants it. The worker is the sole mutator of the authoritative model and the
    // only thread that touches the chunk store; the main thread only ever touches its mirror +
    // meshes.
    public class ChunkManager : MonoBehaviour
    {
        [Header("Grid")]
        public int chunkGridSize = 5;
        public float hexRadius = 1f;
        public int nbIterRelaxation = 2;
        public int nbIterBorderRelaxation = 4;
        public int borderRelaxInteriorRings = 1;
        public bool normalizedRelaxation = false;
        public int seed = 0;

        [Header("Streaming")]
        public int renderRadius = 2;
        [Tooltip("Soft eviction cache: chunks no consumer wants are kept live in memory until " +
                 "the cache exceeds this count, then spilled (FIFO) to the persistent store. " +
                 "A reload while still in the cache costs zero I/O.")]
        public int chunkCacheSize = 32;
        [Tooltip("Warm-presence cache: GameObjects + Meshes + cooked MeshColliders that left the " +
                 "wanted set are kept inactive in the scene until this many accumulate, then " +
                 "FIFO-destroyed. A revisit while still in the cache costs zero BuildMesh / " +
                 "ColliderCook — the major source of camera-pan hitches.")]
        public int presenceCacheSize = 64;
        public Transform cameraTarget;

        [Header("Rendering")]
        [Tooltip("Material applied to the terrain mesh itself (the non-overlay base ground). " +
                 "URP/Lit; the project's terrain shader expects per-vertex barycentric + " +
                 "corner-colour UVs.")]
        public Material groundMaterial;

        [Tooltip("Drives all paintable tile kinds — palette buttons, overlay materials, " +
                 "placement rules, presenter pluggability. Drop a TileCatalog asset here. " +
                 "Adding a new kind = author a new TileDefinition .asset and append it to " +
                 "the catalog's `definitions` list; no edits to this MonoBehaviour or to " +
                 "ChunkSurface required.")]
        public TileCatalog tileCatalog;

        [Header("RTS Highlights")]
        [Tooltip("Colour blended over the cell under the cursor in RTS mode. Alpha controls " +
                 "the blend strength (0.35 is a comfortable default). Pushed to the terrain " +
                 "shader as the global _HoverColor uniform by GameModeController.")]
        public Color hoverColor = new Color(1f, 0.95f, 0.2f, 0.35f);
        [Tooltip("Material applied to the per-chunk outline overlay (line topology, drawn " +
                 "in RTS mode only). URP/Unlit Color is appropriate. Unassigned → outlines " +
                 "render magenta.")]
        public Material outlineMaterial;

        [Header("Elevation")]
        [Tooltip("Tuning knobs for the pointwise terrain height field. Assigned to the static " +
                 "Elevation.Config at Start so worker threads (chunk gen, dual build) and main-" +
                 "thread queries (NPC ground snap) share the same parameters.")]
        public ElevationConfig elevation = new ElevationConfig();

        [Header("Debug")]
        public bool showPrimalGizmos = false;
        public Color primalGizmoColor = new Color(0.3f, 0.8f, 1f);
        public Color edgeVertexColor = new Color(1f, 0.3f, 0.3f);
        public float edgeVertexSize = 0.08f;
        public bool showDualGizmos = false;
        public Color dualGizmoColor = new Color(1f, 0.8f, 0.2f);

        TerrainStreamer streamer;       // owns the authoritative model on a worker thread
        TerrainModel renderModel;       // main-thread mirror the renderer reads from
        ChunkSurface surface;           // owns per-chunk GameObjects (renderer + collider)
        ChunkConsumer renderConsumer;   // this manager's request to the streamer

        // Shared with other systems (e.g. SimulationManager) so they can register their own
        // streaming consumers and collider needs on the surface. Available after Start().
        public TerrainStreamer Streamer => streamer;
        public ChunkSurface Surface => surface;
        public TileCatalog Catalog => tileCatalog;

        // Read-only handle on the main-thread render mirror's view of a chunk. Lets gameplay
        // code (e.g. TilePainter) access published-side data — centroids for nearest-tile
        // lookup, polygon count for bounds — without touching the authoritative worker model.
        // Returns false when the chunk isn't currently in the render set.
        public bool TryGetPublishedChunk(ChunkCoord coord, out PrimalChunk chunk)
            => renderModel.TryGet(coord, out chunk);

        // Exposed scalars so external scripts (TilePainter, save-load) can match the same
        // ChunkCoord math the model uses without duplicating constants.
        public float HexRadius => hexRadius;
        public int ChunkGridSize => chunkGridSize;
        ChunkCoord lastCameraChunk;
        HashSet<ChunkCoord> currentRenderSet; // chunks we want meshed (renderRadius around camera)

        ChunkCoord WorldToChunk(Vector3 pos) => ChunkCoord.FromWorldPos(pos, hexRadius, chunkGridSize);

        HashSet<ChunkCoord> RenderSetAround(ChunkCoord center) =>
            new HashSet<ChunkCoord>(center.HexesInRange(renderRadius));

        void Start()
        {
            if (cameraTarget == null) cameraTarget = Camera.main?.transform;

            // Publish the elevation config before any chunk is generated. Elevation.Sample is
            // pure and reads Config without locking; setting it here once means worker threads
            // see a stable snapshot for the lifetime of the scene.
            Elevation.Config = elevation;

            // Authoritative model — lives behind the streamer, mutated only on the worker.
            TerrainModel model = new TerrainModel
            {
                chunkGridSize = chunkGridSize,
                hexRadius = hexRadius,
                seed = seed,
                nbIterRelaxation = nbIterRelaxation,
                nbIterBorderRelaxation = nbIterBorderRelaxation,
                borderRelaxInteriorRings = borderRelaxInteriorRings,
                normalizedRelaxation = normalizedRelaxation,
                CacheThreshold = chunkCacheSize,
                Store = new MemoryChunkStore(),
            };

            // Mirror — populated with ready chunk copies on the main thread, read by the surface.
            // No store and no generation: the worker is the sole owner of those.
            renderModel = new TerrainModel { chunkGridSize = chunkGridSize, hexRadius = hexRadius };

            // The surface owns the per-chunk GameObjects. We register our render mirror as a
            // primal source; other clients (SimulationManager) add their own mirrors.
            surface = new ChunkSurface(transform, groundMaterial, outlineMaterial,
                                       tileCatalog, hexRadius, chunkGridSize, presenceCacheSize);
            surface.AddPrimalSource(coord => renderModel.TryGet(coord, out PrimalChunk pc) ? pc : null);

            // Register as a render-ready consumer: the streamer loads our requested chunks plus a
            // 1-ring halo, relaxes their seams, and publishes deep copies for us to mesh.
            streamer = new TerrainStreamer(model);
            renderConsumer = streamer.RegisterConsumer(renderReady: true);

            lastCameraChunk = WorldToChunk(cameraTarget.position);
            currentRenderSet = RenderSetAround(lastCameraChunk);
            streamer.SetDesired(renderConsumer, currentRenderSet);

            // Push to the surface with a minimal halo of lower-coord neighbours, so the
            // perimeter cells of currentRenderSet are owned by chunks we actually present
            // (otherwise we'd see seam gaps along the render edge).
            HashSet<ChunkCoord> presented = ChunkCoord.MinimalHalo(currentRenderSet);
            surface.SetRenderRequest(this, presented);
            surface.SetColliderRequest(this, presented);

            // Block once so the opening view isn't empty, then go fully async.
            streamer.Prime();
            DrainResults();
            surface.Apply(); // build the opening meshes immediately so frame 1 isn't blank
            streamer.Start();
        }

        void Update()
        {
            ChunkCoord current = WorldToChunk(cameraTarget.position);
            if (current != lastCameraChunk)
            {
                lastCameraChunk = current;
                currentRenderSet = RenderSetAround(current);
                streamer.SetDesired(renderConsumer, currentRenderSet);

                // The render set is also the collider set (visible terrain must be physical).
                // Both go to the surface expanded by the minimal halo of lower-coord neighbours
                // — they own the perimeter cells under deterministic dual ownership.
                HashSet<ChunkCoord> presented = ChunkCoord.MinimalHalo(currentRenderSet);
                surface.SetRenderRequest(this, presented);
                surface.SetColliderRequest(this, presented);
            }

            DrainResults();
        }

        // The surface reconciles in LateUpdate so it picks up request changes from any other
        // MonoBehaviour (e.g. SimulationManager) that ran earlier or later in the same frame.
        void LateUpdate()
        {
            surface?.Apply();
        }

        // Apply every published pass that arrived since last frame to the mirror model.
        // The surface is reconciled separately in LateUpdate (no need to call it here — a
        // newly-installed chunk will be picked up on the next Apply, ~one frame later).
        void DrainResults()
        {
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.DrainResults);
            while (streamer.TryDrainResult(renderConsumer, out StreamResult result))
            {
                if (result.Error != null) { Debug.LogException(result.Error); continue; }

                using (TerrainProfiler.Measure(TerrainProfiler.Phase.InstallChunk))
                    foreach (PrimalChunk pc in result.Changed)
                        renderModel.Install(pc);

                foreach (ChunkCoord coord in renderModel.LoadedCoords.Where(c => !result.Loaded.Contains(c)).ToList())
                    renderModel.Evict(coord);
            }
        }

        void OnDestroy() => streamer?.Stop();

        void OnDrawGizmos()
        {
            // Gizmos read the main-thread mirror (the worker's model is off-limits here).
            if (renderModel == null) return;

            if (showPrimalGizmos)
            {
                Gizmos.color = primalGizmoColor;
                foreach (ChunkCoord coord in renderModel.LoadedCoords)
                {
                    if (!renderModel.TryGet(coord, out PrimalChunk primal)) continue;

                    // Published render-mirror copies omit the primal graph (DeepCopy nulls
                    // Polygons / Verts — saves the dominant per-publish allocation cost).
                    // Reconstruct primal cell outlines from the CSR neighbour-table arrays
                    // instead: PrimalEdgeVertexA[j] → PrimalEdgeVertexB[j] is one primal
                    // edge, and walking j over each cell's slot range
                    // [PrimalEdgeOffsets[i] .. PrimalEdgeOffsets[i+1]) covers all of that
                    // cell's boundary. Flat chunks ship no neighbour table — silently skip.
                    Vector3[] vAs = primal.PrimalEdgeVertexA;
                    Vector3[] vBs = primal.PrimalEdgeVertexB;
                    if (vAs == null || vBs == null) continue;
                    int n = Mathf.Min(vAs.Length, vBs.Length);
                    for (int j = 0; j < n; j++) Gizmos.DrawLine(vAs[j], vBs[j]);
                }
                // edgeVertexColor / edgeVertexSize spheres were on per-Vertex IsEdge flags;
                // those live on the primal graph that the mirror doesn't ship. Dropped here
                // rather than re-shipping a parallel "is this primal vertex on the chunk
                // border" array just for a gizmo — the cell-outline lines are the part
                // worth seeing.
            }

            if (showDualGizmos && surface != null)
            {
                Gizmos.color = dualGizmoColor;
                foreach (List<Polygon> dual in surface.DualPolygons)
                    foreach (Polygon p in dual) DrawPolygon(p);
            }
        }

        static void DrawPolygon(Polygon p)
        {
            Vector3[] verts = p.GetVerticesPosition();
            for (int i = 0; i < verts.Length; i++)
                Gizmos.DrawLine(verts[i], verts[(i + 1) % verts.Length]);
        }
    }
}
