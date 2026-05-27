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
        public Transform cameraTarget;

        [Header("Rendering")]
        public Material groundMaterial;

        [Header("Debug")]
        public bool showPrimalGizmos = false;
        public Color primalGizmoColor = new Color(0.3f, 0.8f, 1f);
        public Color edgeVertexColor = new Color(1f, 0.3f, 0.3f);
        public float edgeVertexSize = 0.08f;
        public bool showDualGizmos = false;
        public Color dualGizmoColor = new Color(1f, 0.8f, 0.2f);

        TerrainStreamer streamer;       // owns the authoritative model on a worker thread
        TerrainModel renderModel;       // main-thread mirror the renderer reads from
        ChunkRenderer chunkRenderer;
        ChunkConsumer renderConsumer;   // this manager's request to the streamer

        // Shared with other systems (e.g. SimulationManager) so they can register their own
        // consumers. Available after Start(). Null until then.
        public TerrainStreamer Streamer => streamer;
        ChunkCoord lastCameraChunk;
        HashSet<ChunkCoord> currentRenderSet; // chunks we want meshed (renderRadius around camera)

        ChunkCoord WorldToChunk(Vector3 pos) => ChunkCoord.FromWorldPos(pos, hexRadius, chunkGridSize);

        HashSet<ChunkCoord> RenderSetAround(ChunkCoord center) =>
            new HashSet<ChunkCoord>(center.HexesInRange(renderRadius));

        void Start()
        {
            if (cameraTarget == null) cameraTarget = Camera.main?.transform;

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
                Store = new MemoryChunkStore(),
            };

            // Mirror — populated with ready chunk copies on the main thread, read by the renderer.
            // No store and no generation: the worker is the sole owner of those.
            renderModel = new TerrainModel { chunkGridSize = chunkGridSize, hexRadius = hexRadius };
            chunkRenderer = new ChunkRenderer(renderModel, transform, groundMaterial);

            // Register as a render-ready consumer: the streamer loads our requested chunks plus a
            // 1-ring halo, relaxes their seams, and publishes deep copies for us to mesh.
            streamer = new TerrainStreamer(model);
            renderConsumer = streamer.RegisterConsumer(renderReady: true);

            lastCameraChunk = WorldToChunk(cameraTarget.position);
            currentRenderSet = RenderSetAround(lastCameraChunk);
            streamer.SetDesired(renderConsumer, currentRenderSet);
            // Block once so the opening view isn't empty, then go fully async.
            streamer.Prime();
            DrainResults();
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
            }

            DrainResults();
        }

        // Apply every published pass that arrived since last frame to the mirror model, then
        // refresh the render layer once against the current render set.
        void DrainResults()
        {
            bool applied = false;
            while (streamer.TryDrainResult(renderConsumer, out StreamResult result))
            {
                if (result.Error != null) { Debug.LogException(result.Error); continue; }

                foreach (PrimalChunk pc in result.Changed)
                    renderModel.Install(pc);

                foreach (ChunkCoord coord in renderModel.LoadedCoords.Where(c => !result.Loaded.Contains(c)).ToList())
                    renderModel.Evict(coord);

                applied = true;
            }

            if (applied)
                chunkRenderer.SyncTo(currentRenderSet);
        }

        void OnDestroy() => streamer?.Stop();

        void OnDrawGizmos()
        {
            // Gizmos read the main-thread mirror (the worker's model is off-limits here).
            if (renderModel == null) return;

            if (showPrimalGizmos)
            {
                foreach (ChunkCoord coord in renderModel.LoadedCoords)
                {
                    if (!renderModel.TryGet(coord, out PrimalChunk primal)) continue;

                    Gizmos.color = primalGizmoColor;
                    foreach (Polygon p in primal.Polygons) DrawPolygon(p);

                    Gizmos.color = edgeVertexColor;
                    foreach (Polygon p in primal.Polygons)
                        foreach (Vertex v in p.GetVertices())
                            if (v.IsEdge) Gizmos.DrawSphere(v.Position, edgeVertexSize);
                }
            }

            if (showDualGizmos && chunkRenderer != null)
            {
                Gizmos.color = dualGizmoColor;
                foreach (List<Polygon> dual in chunkRenderer.DualPolygons)
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
