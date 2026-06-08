using System;
using System.Collections.Generic;
using TerrainGrid;
using UnityEngine;
using World;

namespace Agents.Economy
{
    // Authoritative map of role-bearing tiles in the world. Survives chunk
    // eviction — a Farmer's AssignedField stays valid even when its chunk
    // is far off-camera. Subscribes to ChunkManager.OnTileChanged so tile
    // mutations propagate into entry creation / role change / destruction;
    // ChunkSurface fires (Default → kind) on chunk first-sight so this
    // registry doesn't need a separate bootstrap scan.
    //
    // Singleton MonoBehaviour: one in the scene, sits on the same
    // "Managers" GameObject as ChunkManager. Awake-time registration is
    // load-bearing — we must subscribe BEFORE ChunkManager.Start runs the
    // initial surface.Apply that emits first-sight events. Unity's
    // guarantee that all Awakes run before any Start makes this safe.
    public class BuildingRegistry : MonoBehaviour
    {
        public static BuildingRegistry Instance { get; private set; }

        [Tooltip("ChunkManager whose tile-change stream feeds the registry. " +
                 "Auto-resolved via FindObjectOfType if left empty.")]
        public ChunkManager chunkManager;

        // Event surface for agent-layer consumers. Idempotent semantics —
        // subscribing later doesn't replay history; query the registry
        // directly for current state.
        public event Action<BuildingEntry> OnEntryAdded;
        public event Action<BuildingEntry> OnEntryRemoved;

        readonly Dictionary<(ChunkCoord, int), BuildingEntry> entries
            = new Dictionary<(ChunkCoord, int), BuildingEntry>();

        // Index by role for cheap "all Markets" / "all Fields" enumeration.
        // Lists are mutated in lockstep with `entries`; no separate
        // invalidation step.
        readonly Dictionary<BuildingRole, List<BuildingEntry>> byRole
            = new Dictionary<BuildingRole, List<BuildingEntry>>();

        void Awake()
        {
            if (Instance != null && Instance != this) { Destroy(this); return; }
            Instance = this;

            if (chunkManager == null) chunkManager = FindObjectOfType<ChunkManager>();
            if (chunkManager == null)
            {
                Debug.LogError("BuildingRegistry: no ChunkManager in scene; disabling.");
                enabled = false;
                return;
            }

            // Subscribe in Awake so chunks that load during ChunkManager.Start's
            // initial surface.Apply find us already listening.
            chunkManager.OnTileChanged += HandleTileChanged;
        }

        void OnEnable()
        {
            // Defer the clock subscription to OnEnable because the GameClock
            // instance might not exist yet at Awake time (different script
            // execution orders are tolerated this way).
            if (GameClock.Instance != null) GameClock.Instance.OnNewDay += HandleNewDay;
        }

        void OnDisable()
        {
            if (GameClock.Instance != null) GameClock.Instance.OnNewDay -= HandleNewDay;
        }

        void OnDestroy()
        {
            if (chunkManager != null) chunkManager.OnTileChanged -= HandleTileChanged;
            if (Instance == this) Instance = null;
        }

        // -------- Tile-change handling --------

        // Idempotent: the same (coord, cell, kind) being reported repeatedly
        // (e.g. a chunk parks, gets evicted, reloads → first-sight events
        // replay) is a no-op as long as the *current registry role* matches
        // the new event's role. The event's `oldK` is advisory and ignored —
        // we trust our own state.
        void HandleTileChanged(ChunkCoord coord, int cell, TileKind oldK, TileKind newK)
        {
            BuildingRole newRole = BuildingRoles.FromTileKind(newK);
            BuildingRole curRole = entries.TryGetValue((coord, cell), out BuildingEntry existing)
                ? existing.Role : BuildingRole.None;

            if (curRole == newRole) return;

            if (curRole != BuildingRole.None) RemoveEntryInternal(coord, cell);
            if (newRole != BuildingRole.None) AddEntryInternal(coord, cell, newRole);
        }

        void AddEntryInternal(ChunkCoord coord, int cell, BuildingRole role)
        {
            Vector3 worldPos = LookupCentroid(coord, cell);
            int slots = BuildingRoles.DefaultSlots(role);
            BuildingEntry e = new BuildingEntry(coord, cell, role, worldPos, slots);
            if (role == BuildingRole.Market) e.Market = new MarketStockpile();

            entries[(coord, cell)] = e;
            if (!byRole.TryGetValue(role, out List<BuildingEntry> list))
                byRole[role] = list = new List<BuildingEntry>();
            list.Add(e);
            OnEntryAdded?.Invoke(e);
        }

        void RemoveEntryInternal(ChunkCoord coord, int cell)
        {
            if (!entries.TryGetValue((coord, cell), out BuildingEntry e)) return;
            entries.Remove((coord, cell));
            if (byRole.TryGetValue(e.Role, out List<BuildingEntry> list)) list.Remove(e);
            // Release any squatter claims so dependent agents replan instead
            // of holding ghost references.
            e.SlotsClaimed.Clear();
            OnEntryRemoved?.Invoke(e);
        }

        // Centroid lookup against the published mirror. Returns Vector3.zero
        // if the chunk isn't currently loaded — agents that target this
        // entry will need to either wait or replan. The position update on
        // chunk reload would be a nice future hardening; for v1 the chunk
        // is loaded when the player paints, so this is always reachable.
        Vector3 LookupCentroid(ChunkCoord coord, int cell)
        {
            if (!chunkManager.TryGetPublishedChunk(coord, out PrimalChunk pc)) return Vector3.zero;
            if (pc.PrimalCellCentroids == null || cell < 0 || cell >= pc.PrimalCellCentroids.Length)
                return Vector3.zero;
            return pc.PrimalCellCentroids[cell];
        }

        // -------- Clock-driven price tick --------

        void HandleNewDay(int day)
        {
            if (!byRole.TryGetValue(BuildingRole.Market, out List<BuildingEntry> markets)) return;
            for (int i = 0; i < markets.Count; i++)
            {
                MarketStockpile m = markets[i].Market;
                if (m != null) PriceModel.DailyTick(m);
            }
        }

        // -------- Public query API --------

        public bool TryGet(ChunkCoord coord, int cell, out BuildingEntry entry)
            => entries.TryGetValue((coord, cell), out entry);

        public IReadOnlyList<BuildingEntry> AllOfRole(BuildingRole role)
        {
            if (byRole.TryGetValue(role, out List<BuildingEntry> list)) return list;
            return System.Array.Empty<BuildingEntry>();
        }

        public int CountOfRole(BuildingRole role)
            => byRole.TryGetValue(role, out List<BuildingEntry> list) ? list.Count : 0;

        // Brute-force nearest scan. Cheap at v1 sizes (≤ a few hundred role-
        // bearing tiles) and avoids the index-invalidation hazards of a
        // spatial hash that the entries would need to be updated in. Add a
        // spatial index later if profile shows need.
        public BuildingEntry FindNearestByRole(Vector3 pos, BuildingRole role)
        {
            if (!byRole.TryGetValue(role, out List<BuildingEntry> list) || list.Count == 0) return null;
            BuildingEntry best = null;
            float bestSq = float.PositiveInfinity;
            for (int i = 0; i < list.Count; i++)
            {
                float d = (list[i].WorldPos - pos).sqrMagnitude;
                if (d < bestSq) { bestSq = d; best = list[i]; }
            }
            return best;
        }

        // -------- Slot reservation --------

        // First-come-first-served claim. Owners are opaque to the registry —
        // pass `this` from an Agent component, or any object that identifies
        // the claimant. Returns true on success; on failure the caller
        // should treat the action as IsValid()==false and replan.
        public bool TryClaimSlot(BuildingEntry e, object owner)
        {
            if (e == null || owner == null) return false;
            if (e.SlotsClaimed.Contains(owner)) return true;   // re-claim is a no-op
            if (!e.HasFreeSlot) return false;
            e.SlotsClaimed.Add(owner);
            return true;
        }

        public void ReleaseSlot(BuildingEntry e, object owner)
        {
            if (e == null || owner == null) return;
            e.SlotsClaimed.Remove(owner);
        }

        // Release every slot this owner is holding — useful when an Agent
        // dies / despawns mid-plan and we can't trust it to call Release
        // for itself.
        public void ReleaseAllForOwner(object owner)
        {
            if (owner == null) return;
            foreach (BuildingEntry e in entries.Values)
                e.SlotsClaimed.Remove(owner);
        }

        // -------- Debug HUD --------

        [Header("Debug")]
        public bool drawHud = false;
        public KeyCode toggleHudKey = KeyCode.F4;

        void Update()
        {
            if (Input.GetKeyDown(toggleHudKey)) drawHud = !drawHud;
        }

        void OnGUI()
        {
            if (!drawHud) return;
            const float w = 260f;
            float x = Screen.width - w - 10f;
            float y = 70f;
            GUI.Box(new Rect(x, y, w, 22f * (1 + 5)), "");
            GUI.Label(new Rect(x + 6f, y + 2f, w - 12f, 20f), $"BuildingRegistry  [{toggleHudKey} hide]");
            int row = 1;
            foreach (BuildingRole r in (BuildingRole[])System.Enum.GetValues(typeof(BuildingRole)))
            {
                if (r == BuildingRole.None) continue;
                int n = CountOfRole(r);
                GUI.Label(new Rect(x + 6f, y + 2f + row * 18f, w - 12f, 18f),
                          $"  {BuildingRoles.Name(r)}: {n}");
                row++;
            }
        }
    }
}
