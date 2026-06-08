using System.Collections.Generic;
using Agents.Core;
using Agents.Economy;
using Agents.Pathfinding;
using TerrainGrid;
using UnityEngine;

namespace Agents.Runtime
{
    // Demo bootstrap: paint a village (1+ House, 1+ Field, 1+ Bakery,
    // 1+ Market) and press the spawn key. The spawner reads the
    // BuildingRegistry, distributes Farmers / Bakers / Citizens across
    // the available tiles, and instantiates each as an Agent.
    //
    // Why manual key trigger instead of auto-spawn: gives the player
    // control over WHEN the demo starts. Some scenes need a few seconds
    // of camera setup or terrain streaming before a paint pass; auto-
    // spawn at Start would race those.
    //
    // Limitations (v1):
    //   • No de-duplication — each press of the spawn key adds another
    //     batch. Useful for population growth tests; reset the scene
    //     between runs if you want a clean count.
    //   • Citizens with no Job tile get a starter cash float so they
    //     can sustain RestockGoal until they run out — they DO starve
    //     long-term. Future: passive income source / employment.
    public class VillageSpawner : MonoBehaviour
    {
        [Header("Refs")]
        [Tooltip("Auto-found via FindObjectOfType if left empty.")]
        public ChunkManager chunkManager;
        public GameObject agentPrefab;

        [Header("Roster")]
        public int farmers = 2;
        public int bakers  = 1;
        public int citizens = 2;

        [Tooltip("Starting cash for Citizens — they pay for bread until they run out.")]
        public int citizenStartMoney = 200;

        [Tooltip("Starting cash for Bakers — they need to buy Wheat at the Market before they " +
                 "can produce Bread. Without this float their EarnGoal can't form a plan and " +
                 "they idle.")]
        public int bakerStartMoney = 50;

        [Tooltip("Starting cash for Farmers. Zero is fine — they earn from selling Wheat.")]
        public int farmerStartMoney = 0;

        [Tooltip("Press this to instantiate the roster. The BuildingRegistry must already " +
                 "contain at least one Market and one House.")]
        public KeyCode spawnKey = KeyCode.F9;

        [Tooltip("Auto-spawn at Start once the registry has the required tiles. Off by default " +
                 "because the player typically needs a beat to paint the village first.")]
        public bool autoSpawnIfReady = false;

        GoapBlueprint blueprint;
        TilePathCache pathCache;
        int spawnedSoFar;

        void Start()
        {
            if (chunkManager == null) chunkManager = FindObjectOfType<ChunkManager>();
            if (chunkManager == null)
            {
                Debug.LogError("VillageSpawner: no ChunkManager; disabling.");
                enabled = false;
                return;
            }
            pathCache = new TilePathCache(chunkManager, 128);
            blueprint = VillageBlueprint.Build(chunkManager, pathCache);

            if (autoSpawnIfReady && IsReady()) Spawn();
        }

        void Update()
        {
            if (Input.GetKeyDown(spawnKey)) Spawn();
        }

        bool IsReady()
        {
            BuildingRegistry reg = BuildingRegistry.Instance;
            return reg != null
                && reg.CountOfRole(BuildingRole.House)  >= 1
                && reg.CountOfRole(BuildingRole.Market) >= 1;
        }

        public void Spawn()
        {
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) { Debug.LogWarning("VillageSpawner: no BuildingRegistry."); return; }

            IReadOnlyList<BuildingEntry> houses   = reg.AllOfRole(BuildingRole.House);
            IReadOnlyList<BuildingEntry> fields   = reg.AllOfRole(BuildingRole.Field);
            IReadOnlyList<BuildingEntry> bakeries = reg.AllOfRole(BuildingRole.Bakery);
            IReadOnlyList<BuildingEntry> markets  = reg.AllOfRole(BuildingRole.Market);

            if (markets.Count == 0)
            {
                Debug.LogWarning("VillageSpawner: paint at least one Market before spawning.");
                return;
            }
            if (houses.Count == 0)
            {
                Debug.LogWarning("VillageSpawner: paint at least one House before spawning.");
                return;
            }

            int idx = 0;

            // Farmers — only spawned if there are Fields to assign.
            int farmerCount = fields.Count > 0 ? farmers : 0;
            for (int i = 0; i < farmerCount; i++)
                SpawnOne("Farmer", houses[idx % houses.Count], fields[i % fields.Count], markets[0], farmerStartMoney, ref idx);

            // Bakers — only spawned if there are Bakeries.
            int bakerCount = bakeries.Count > 0 ? bakers : 0;
            for (int i = 0; i < bakerCount; i++)
                SpawnOne("Baker", houses[idx % houses.Count], bakeries[i % bakeries.Count], markets[0], bakerStartMoney, ref idx);

            // Citizens — no Job assignment.
            for (int i = 0; i < citizens; i++)
                SpawnOne("Citizen", houses[idx % houses.Count], null, markets[0], citizenStartMoney, ref idx);
        }

        void SpawnOne(string roleName, BuildingEntry home, BuildingEntry job,
                      BuildingEntry market, int startMoney, ref int idx)
        {
            Vector3 pos = home.WorldPos + new Vector3(0f, 0.5f, 0f);
            GameObject go = agentPrefab != null
                ? Instantiate(agentPrefab, pos, Quaternion.identity)
                : CreateFallbackVisual(pos, roleName, idx);

            Agent agent = go.GetComponent<Agent>();
            if (agent == null) agent = go.AddComponent<Agent>();

            (ChunkCoord, int)? jobTuple = null;
            if (job != null) jobTuple = (job.Coord, job.CellIndex);

            agent.Initialize(blueprint,
                             $"{roleName}{spawnedSoFar}",
                             (home.Coord, home.CellIndex),
                             jobTuple);

            // Seed initial target keys in symbolic state so sensors /
            // actions can resolve assignments without a separate
            // bootstrap pass. Stored as plain value tuples — sensors
            // retrieve them via Get<(ChunkCoord, int)>().
            agent.State.Set(TargetKey.Home, (home.Coord, home.CellIndex));
            if (job != null)
            {
                if (job.Role == BuildingRole.Field)
                    agent.State.Set(TargetKey.AssignedField, (job.Coord, job.CellIndex));
                else if (job.Role == BuildingRole.Bakery)
                    agent.State.Set(TargetKey.AssignedBakery, (job.Coord, job.CellIndex));
            }
            agent.State.Set(TargetKey.NearestMarket, (market.Coord, market.CellIndex));

            agent.Money = startMoney;
            spawnedSoFar++;
            idx++;
        }

        static GameObject CreateFallbackVisual(Vector3 pos, string roleName, int idx)
        {
            GameObject go = GameObject.CreatePrimitive(PrimitiveType.Capsule);
            go.name = $"{roleName}{idx}";
            go.transform.position = pos;
            Collider col = go.GetComponent<Collider>();
            if (col != null) Destroy(col);
            // Tint by role for at-a-glance identification.
            MeshRenderer mr = go.GetComponent<MeshRenderer>();
            if (mr != null)
            {
                Color tint = roleName switch
                {
                    "Farmer"  => new Color(0.5f, 0.8f, 0.3f),
                    "Baker"   => new Color(0.9f, 0.7f, 0.3f),
                    "Citizen" => new Color(0.6f, 0.6f, 0.9f),
                    _         => Color.white,
                };
                mr.material.color = tint;
            }
            return go;
        }
    }
}
