using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid
{
    // A first NPC simulator and a second consumer of the terrain streamer (alongside the
    // renderer). It registers a *primal-only* consumer — it wants the cheap data layer around a
    // focus point, never the dual/mesh — and keeps a crowd of agents wandering within that bubble
    // in continuous world space.
    //
    // Scaling intent baked in from the start:
    //   • Agents are struct data in a flat array, updated in batch (no MonoBehaviour per NPC).
    //   • Simulation is decoupled from rendering: every agent is ticked, but only agents near the
    //     camera are drawn (GPU-instanced). Dormant/fast-forward tiers are a later phase.
    //   • The terrain dependency goes through a private simModel mirror, exactly like the
    //     renderer's renderModel — thread-safe, never touching the worker's authoritative model.
    public class SimulationManager : MonoBehaviour
    {
        [Header("Terrain")]
        public ChunkManager terrain;     // source of the streamer + grid config
        public Transform focus;          // streaming focus (defaults to terrain.cameraTarget)

        [Header("Simulation")]
        public int agentCount = 30;
        public int simRadius = 3;        // chunks around focus kept loaded (primal-only)
        public float agentSpeed = 3f;
        public float retargetDistance = 1f;
        public int seed = 12345;

        [Header("Rendering")]
        public Mesh agentMesh;           // e.g. a capsule/cube; material must enable GPU instancing
        public Material agentMaterial;
        public Vector3 agentScale = new Vector3(0.5f, 1f, 0.5f);
        public float agentYOffset = 1f;  // lift the mesh so its base sits on the ground
        public float renderRadiusWorld = 60f;

        TerrainModel simModel;           // main-thread mirror populated from the sim consumer
        ChunkConsumer consumer;
        Agent[] agents;
        System.Random rng;
        ChunkCoord lastFocusChunk;
        bool hasFocusChunk;

        Matrix4x4[] instanceBuffer;
        bool warnedMissingVisual;
        bool spawned;

        float ChunkSize => terrain.chunkGridSize * terrain.hexRadius;
        // Soft world radius of the sim bubble (flat-top hex circumradius times the ring count).
        float BubbleRadiusWorld => (simRadius + 0.5f) * ChunkSize;

        void Start()
        {
            if (terrain == null) terrain = FindFirstObjectByType<ChunkManager>();
            if (terrain == null)
            {
                Debug.LogError("SimulationManager: no ChunkManager found; disabling.");
                enabled = false;
                return;
            }
            simModel = new TerrainModel { chunkGridSize = terrain.chunkGridSize, hexRadius = terrain.hexRadius };
            rng = new System.Random(seed);
            instanceBuffer = new Matrix4x4[1023]; // GPU instancing draws at most 1023 per call

            // RenderMeshInstanced requires the material to opt into GPU instancing; enable it
            // here so it works regardless of the asset's checkbox.
            if (agentMaterial != null) agentMaterial.enableInstancing = true;
        }

        // Register lazily: the streamer is created in ChunkManager.Start, and Unity doesn't
        // guarantee Start order, so we wait until it exists.
        void EnsureConsumer()
        {
            if (consumer != null) return;
            TerrainStreamer streamer = terrain != null ? terrain.Streamer : null;
            if (streamer != null) consumer = streamer.RegisterConsumer(renderReady: false);
        }

        void Update()
        {
            if (focus == null && terrain != null) focus = terrain.cameraTarget;
            EnsureConsumer();
            if (consumer == null || focus == null) return;

            if (!spawned) { SpawnAgents(); spawned = true; }

            UpdateDesired();
            DrainSimResults();
            TickAgents(Time.deltaTime);
            RenderAgents();
        }

        void OnDestroy()
        {
            if (consumer != null) terrain?.Streamer?.UnregisterConsumer(consumer);
        }

        void SpawnAgents()
        {
            agents = new Agent[Mathf.Max(0, agentCount)];
            Vector3 center = focus != null ? Flat(focus.position) : Vector3.zero;
            for (int i = 0; i < agents.Length; i++)
            {
                Vector3 p = RandomPointInBubble(center);
                agents[i] = new Agent
                {
                    Position = Ground(p),
                    Velocity = Vector3.zero,
                    Target = RandomPointInBubble(center),
                    Speed = agentSpeed,
                };
            }
        }

        // Declare the primal-only chunks we want loaded around the focus (only when it changes).
        void UpdateDesired()
        {
            ChunkCoord c = ChunkCoord.FromWorldPos(focus.position, terrain.hexRadius, terrain.chunkGridSize);
            if (hasFocusChunk && c == lastFocusChunk) return;
            lastFocusChunk = c;
            hasFocusChunk = true;
            terrain.Streamer.SetDesired(consumer, new HashSet<ChunkCoord>(c.HexesInRange(simRadius)));
        }

        // Mirror the published primal copies into simModel (same pattern as ChunkManager).
        void DrainSimResults()
        {
            while (terrain.Streamer.TryDrainResult(consumer, out StreamResult result))
            {
                if (result.Error != null) { Debug.LogException(result.Error); continue; }

                foreach (PrimalChunk pc in result.Changed)
                    simModel.Install(pc);

                foreach (ChunkCoord coord in new List<ChunkCoord>(simModel.LoadedCoords))
                    if (!result.Loaded.Contains(coord)) simModel.Evict(coord);
            }
        }

        void TickAgents(float dt)
        {
            Vector3 center = Flat(focus.position);
            float bound = BubbleRadiusWorld;

            for (int i = 0; i < agents.Length; i++)
            {
                ref Agent a = ref agents[i];

                Vector3 toTarget = Flat(a.Target - a.Position);
                if (toTarget.magnitude < retargetDistance)
                {
                    a.Target = RandomPointInBubble(center);
                    toTarget = Flat(a.Target - a.Position);
                }

                // Simple steering toward the target.
                Vector3 desired = toTarget.sqrMagnitude > 1e-4f ? toTarget.normalized * a.Speed : Vector3.zero;
                a.Velocity = Vector3.MoveTowards(a.Velocity, desired, a.Speed * 4f * dt);
                a.Position += a.Velocity * dt;

                // Keep inside the loaded bubble (it tracks the camera, so agents stay on terrain).
                Vector3 fromCenter = Flat(a.Position - center);
                if (fromCenter.magnitude > bound)
                {
                    a.Position = center + fromCenter.normalized * bound;
                    a.Target = RandomPointInBubble(center);
                }

                a.Position = Ground(a.Position);
            }
        }

        void RenderAgents()
        {
            if (agentMesh == null || agentMaterial == null)
            {
                if (!warnedMissingVisual)
                {
                    Debug.LogWarning("SimulationManager: assign agentMesh and a GPU-instancing material to see agents.");
                    warnedMissingVisual = true;
                }
                return;
            }

            Vector3 camPos = Flat(focus.position);
            float renderSqr = renderRadiusWorld * renderRadiusWorld;
            var rp = new RenderParams(agentMaterial);

            int count = 0;
            for (int i = 0; i < agents.Length; i++)
            {
                if (Flat(agents[i].Position - camPos).sqrMagnitude > renderSqr) continue; // tier: cull far agents

                Vector3 vel = Flat(agents[i].Velocity);
                Quaternion rot = vel.sqrMagnitude > 1e-4f ? Quaternion.LookRotation(vel.normalized, Vector3.up) : Quaternion.identity;
                instanceBuffer[count++] = Matrix4x4.TRS(agents[i].Position, rot, agentScale);

                if (count == instanceBuffer.Length)
                {
                    Graphics.RenderMeshInstanced(rp, agentMesh, 0, instanceBuffer, count);
                    count = 0;
                }
            }
            if (count > 0) Graphics.RenderMeshInstanced(rp, agentMesh, 0, instanceBuffer, count);
        }

        // --- helpers ---

        static Vector3 Flat(Vector3 v) => new Vector3(v.x, 0f, v.z);

        Vector3 RandomPointInBubble(Vector3 center)
        {
            float ang = (float)rng.NextDouble() * Mathf.PI * 2f;
            float r = Mathf.Sqrt((float)rng.NextDouble()) * BubbleRadiusWorld; // uniform over disc
            return center + new Vector3(Mathf.Cos(ang) * r, 0f, Mathf.Sin(ang) * r);
        }

        // Snap a world position onto the ground. Flat (Y=0) for now; this is the hook that will
        // sample the primal height field / terrain surface once it's non-flat.
        Vector3 Ground(Vector3 p) => new Vector3(p.x, GroundHeight(p) + agentYOffset, p.z);
        float GroundHeight(Vector3 p) => 0f;
    }
}
