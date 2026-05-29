using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Pool;

namespace TerrainGrid
{
    // The NPC crowd controller. Three responsibilities:
    //   1. As a streaming consumer, keep the chunks under each agent loaded (plus a small
    //      buffer ring around each one for the agent's next step, and a buffer around `focus`
    //      to guarantee somewhere to spawn). The desired set is recomputed when the agents'
    //      chunk membership changes.
    //   2. As a collider client of ChunkManager.Surface, request MeshColliders on those same
    //      chunks so the agents have something to walk on — *without* attaching a renderer
    //      when the chunk is outside the camera's render set.
    //   3. Spawn a pool of physics-driven NPCAgents inside the focus bubble, give each one a
    //      wander target, and recycle any that drift too far away.
    //
    // The agents are real GameObjects with Rigidbody + Collider, so the boat (or anything
    // else) can collide with them and push them around. Object pooling (UnityEngine.Pool)
    // is used from day one: spawn = SetActive(true), despawn = SetActive(false).
    //
    // The streaming consumer is renderReady=true because we need dual-capable chunks to
    // build the collider mesh — primal-only chunks don't have relaxed seams and would crack.
    public class SimulationManager : MonoBehaviour
    {
        [Header("Terrain")]
        public ChunkManager terrain;     // source of the streamer + grid config
        public Transform focus;          // streaming focus (defaults to terrain.cameraTarget)

        [Header("Simulation")]
        public NPCAgent agentPrefab;     // must have Rigidbody + Collider + a renderer
        public int agentCount = 200;
        [Tooltip("Soft bubble radius (in chunks around focus) that bounds where agents wander and respawn.")]
        public int simRadius = 4;
        [Tooltip("Extra rings of chunks loaded around each agent (1 = chunk under agent + 6 neighbours).")]
        public int agentBuffer = 1;
        [Tooltip("Rings of chunks kept loaded around focus regardless of agents, so there's always terrain to spawn onto.")]
        public int spawnBuffer = 1;
        public float agentSpeed = 3f;
        public float retargetDistance = 1f;
        public float spawnHeight = 2f;   // initial Y when spawning (above ground)
        public int seed = 12345;

        [Tooltip("Cap on how many NPCs are instantiated per frame on startup. " +
                 "Instantiating a physics prefab costs ~0.3-1ms each, so spawning the whole " +
                 "crowd at once would visibly hitch. Pool growth during play uses the same cap.")]
        public int spawnsPerFrame = 20;

        TerrainModel simModel;           // main-thread mirror populated from the sim consumer
        ChunkConsumer consumer;
        ChunkSurface surface;            // borrowed from ChunkManager; we add collider requests to it
        System.Random rng;
        HashSet<ChunkCoord> lastDesired; // last set sent to streamer + surface (for change detection)
        HashSet<ChunkCoord> scratch = new HashSet<ChunkCoord>(); // reused per Update to dodge alloc

        ObjectPool<NPCAgent> pool;
        readonly List<NPCAgent> active = new();
        int spawnedSoFar;                 // initial-spawn progress (staggered across frames)
        bool warnedMissingPrefab;

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

            pool = new ObjectPool<NPCAgent>(
                createFunc: () =>
                {
                    NPCAgent a = Instantiate(agentPrefab, transform);
                    a.gameObject.SetActive(false);
                    return a;
                },
                actionOnGet: a => a.gameObject.SetActive(true),
                actionOnRelease: a => a.gameObject.SetActive(false),
                actionOnDestroy: a => { if (a != null) Destroy(a.gameObject); },
                collectionCheck: false,
                // Defaults sized for agentCount, with headroom for churn (a small wave of
                // recycled agents in flight). Pool never destroys instances under maxSize,
                // so steady-state operation is allocation-free.
                defaultCapacity: Mathf.Max(agentCount, 8),
                maxSize: Mathf.Max(agentCount + 32, 64)
            );
        }

        // Wire up lazily: the streamer and surface are created in ChunkManager.Start, and
        // Unity doesn't guarantee Start order, so we wait until they exist.
        void EnsureLinks()
        {
            if (consumer == null && terrain?.Streamer != null)
                consumer = terrain.Streamer.RegisterConsumer(renderReady: true);

            if (surface == null && terrain?.Surface != null)
            {
                surface = terrain.Surface;
                // Register our mirror as a primal source so the surface can build collider
                // meshes for chunks the camera doesn't render.
                surface.AddPrimalSource(coord => simModel.TryGet(coord, out PrimalChunk pc) ? pc : null);
            }
        }

        void Update()
        {
            if (focus == null && terrain != null) focus = terrain.cameraTarget;
            EnsureLinks();
            if (consumer == null || focus == null) return;

            if (agentPrefab == null)
            {
                if (!warnedMissingPrefab)
                {
                    Debug.LogWarning("SimulationManager: assign an agentPrefab (with Rigidbody + Collider) to see NPCs.");
                    warnedMissingPrefab = true;
                }
                return;
            }

            UpdateDesired();       // recompute agent-driven chunk set, push to streamer + surface
            DrainSimResults();
            SpawnInitialBatch();   // staggered — runs across multiple frames until full
            UpdateAgents();
        }

        void OnDestroy()
        {
            if (consumer != null) terrain?.Streamer?.UnregisterConsumer(consumer);
            surface?.SetColliderRequest(this, null); // drop our collider claim
            pool?.Dispose();
        }

        // Drip-spawn the initial crowd over multiple frames so we don't pay 200×Instantiate
        // in a single frame (which would visibly hitch). Once spawnedSoFar reaches agentCount,
        // this is a no-op for the rest of the session.
        void SpawnInitialBatch()
        {
            if (spawnedSoFar >= agentCount) return;
            Vector3 center = Flat(focus.position);
            int budget = Mathf.Min(spawnsPerFrame, agentCount - spawnedSoFar);
            for (int i = 0; i < budget; i++) SpawnOne(center);
            spawnedSoFar += budget;
        }

        void SpawnOne(Vector3 center)
        {
            NPCAgent a = pool.Get();
            Vector3 p = RandomPointInBubble(center);
            a.transform.SetPositionAndRotation(
                new Vector3(p.x, GroundHeight(p) + spawnHeight, p.z),
                Quaternion.identity);

            // Reset any leftover physics state from a previous pool use.
            if (a.TryGetComponent(out Rigidbody rb))
            {
                rb.linearVelocity = Vector3.zero;
                rb.angularVelocity = Vector3.zero;
            }

            a.speed = agentSpeed;
            a.Target = RandomPointInBubble(center);
            active.Add(a);
        }

        // The chunks we want loaded = (chunks under each agent, expanded by agentBuffer rings)
        // ∪ (chunks around focus, expanded by spawnBuffer rings). Recomputed every frame; we
        // only re-publish to the streamer + surface if the set actually changed, so static
        // crowds cause no churn.
        void UpdateDesired()
        {
            scratch.Clear();

            // Always keep a small ring around focus loaded so there's terrain to spawn onto.
            ChunkCoord focusCoord = ChunkCoord.FromWorldPos(focus.position, terrain.hexRadius, terrain.chunkGridSize);
            foreach (ChunkCoord c in focusCoord.HexesInRange(spawnBuffer))
                scratch.Add(c);

            // One chunk per agent + buffer ring (so the agent's next step has ground waiting).
            for (int i = 0; i < active.Count; i++)
            {
                Vector3 p = active[i].transform.position;
                ChunkCoord ac = ChunkCoord.FromWorldPos(p, terrain.hexRadius, terrain.chunkGridSize);
                foreach (ChunkCoord c in ac.HexesInRange(agentBuffer))
                    scratch.Add(c);
            }

            if (lastDesired != null && scratch.SetEquals(lastDesired)) return;

            lastDesired = new HashSet<ChunkCoord>(scratch);
            terrain.Streamer.SetDesired(consumer, lastDesired);

            // Push to the surface with a minimal halo: the lower-coord neighbours of every
            // chunk in lastDesired also need a collider, because they own the perimeter cells
            // along the shared seams. Without this we'd see hairline gaps under agents that
            // straddle a chunk edge.
            surface?.SetColliderRequest(this, ChunkCoord.MinimalHalo(lastDesired));
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

        // Crowd-level bookkeeping: pick new targets, recycle stragglers. Steering itself is
        // owned by each NPCAgent's FixedUpdate, so this only writes intent (Target).
        void UpdateAgents()
        {
            Vector3 center = Flat(focus.position);
            float bound = BubbleRadiusWorld;
            float boundSqr = bound * bound;

            for (int i = active.Count - 1; i >= 0; i--)
            {
                NPCAgent a = active[i];
                Vector3 pos = a.transform.position;

                // If a physics object has fallen through the world (or generally gone rogue),
                // pull it back into the pool and respawn fresh.
                if (pos.y < -50f || Flat(pos - center).sqrMagnitude > boundSqr * 4f)
                {
                    active.RemoveAt(i);
                    pool.Release(a);
                    SpawnOne(center);
                    continue;
                }

                if (a.ReachedTarget(retargetDistance))
                    a.Target = RandomPointInBubble(center);

                // Soft containment: if the agent wandered past the loaded bubble, send it back
                // toward the center instead of letting it walk off the terrain.
                if (Flat(pos - center).sqrMagnitude > boundSqr)
                    a.Target = center;
            }
        }

        // --- helpers ---

        static Vector3 Flat(Vector3 v) => new Vector3(v.x, 0f, v.z);

        Vector3 RandomPointInBubble(Vector3 center)
        {
            float ang = (float)rng.NextDouble() * Mathf.PI * 2f;
            float r = Mathf.Sqrt((float)rng.NextDouble()) * BubbleRadiusWorld; // uniform over disc
            return center + new Vector3(Mathf.Cos(ang) * r, 0f, Mathf.Sin(ang) * r);
        }

        // Hook for sampling the primal height field once terrain stops being flat. Until then,
        // ground is Y=0 and the agent's collider + gravity handle vertical placement.
        float GroundHeight(Vector3 p) => 0f;
    }
}
