using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace TerrainGrid
{
    // The abstract, non-rendering terrain: owns every loaded primal chunk and the
    // data-layer operations (generation + cross-chunk border relaxation). Knows nothing
    // about cameras or meshes, so it can be queried for gameplay (occupancy, adjacency)
    // and — being free of Unity rendering — is a natural candidate to move off-thread.
    //
    // The render layer reads this model; this model never references the render layer.
    public class TerrainModel
    {
        public int chunkGridSize = 5;
        public float hexRadius = 1f;
        public int seed = 0;
        public int nbIterRelaxation = 2;
        public int nbIterBorderRelaxation = 4;
        public int borderRelaxInteriorRings = 1;
        public bool normalizedRelaxation = false;

        // Persistent store: unloaded chunks are snapshotted here and restored on reload, so
        // their relaxed state survives instead of being regenerated fresh. May be null
        // (then unload truly discards and reload regenerates).
        public IChunkStore Store;

        // Soft cap on the in-memory eviction cache (chunks no consumer wants but kept live
        // for cheap revival). When the cache exceeds this, the oldest entries are flushed to
        // the store (FIFO). Set to 0 to skip the cache entirely.
        public int CacheThreshold = 32;

        readonly Dictionary<ChunkCoord, PrimalChunk> chunks = new();

        // Shared empty-list sentinel for flat chunks' Dual: lets the version machinery treat
        // them as "dual built (zero polygons)" without allocating a fresh list per chunk.
        static readonly List<Polygon> noPolygons = new List<Polygon>();

        // Soft-eviction tier: chunks that no consumer wants right now but are kept in memory
        // so the next pass can promote them back with no I/O. Promotion happens via
        // TryReactivate; overflow flushes to Store via TrimCache. cacheOrder is FIFO and may
        // contain stale entries (promoted-out coords), which are skipped during trimming.
        readonly Dictionary<ChunkCoord, PrimalChunk> cache = new();
        readonly Queue<ChunkCoord> cacheOrder = new();

        // Tile mutations targeting chunks that weren't loaded when they were enqueued. Replayed
        // by ReplayPendingMutations the moment that coord comes live (generate / cache revive /
        // store restore), so gameplay code never has to wait for streaming to apply a paint.
        readonly Dictionary<ChunkCoord, List<TileMutation>> pendingMutations = new();

        // Guards every read/write of `pendingMutations` and of any chunk's `TileProperties`
        // array. Taken by ApplyTileMutation, ReplayPendingMutations, and TryGetTileProperty —
        // the three paths that touch tile state. Held only for tight, allocation-free spans so
        // main-thread queries never block the worker meaningfully.
        readonly object tileLock = new object();

        public IEnumerable<ChunkCoord> LoadedCoords => chunks.Keys;
        public int LoadedCount => chunks.Count;
        public int CachedCount => cache.Count;
        public bool IsLoaded(ChunkCoord coord) => chunks.ContainsKey(coord);
        public bool IsCached(ChunkCoord coord) => cache.ContainsKey(coord);
        public bool TryGet(ChunkCoord coord, out PrimalChunk chunk) => chunks.TryGetValue(coord, out chunk);

        // --- Tile property mutation (worker-thread) ---
        // All writes go through ApplyTileMutation; main thread calls in via TerrainStreamer.
        // If the target chunk is loaded, write straight to its lazy-allocated parallel array
        // and bump Version (existing publish loop catches the change). If it's not loaded,
        // park the mutation in pendingMutations so the next promotion (generate/revive/restore)
        // applies it via ReplayPendingMutations before any consumer sees the chunk.

        public void ApplyTileMutation(ChunkCoord coord, int polygonIndex, TileProperty value)
        {
            lock (tileLock)
            {
                if (chunks.TryGetValue(coord, out PrimalChunk chunk))
                {
                    if (polygonIndex < 0 || polygonIndex >= chunk.Polygons.Count) return;
                    if (chunk.TileProperties == null)
                        chunk.TileProperties = new TileProperty[chunk.Polygons.Count];
                    chunk.TileProperties[polygonIndex] = value;
                    chunk.Version++;
                    return;
                }

                if (!pendingMutations.TryGetValue(coord, out List<TileMutation> list))
                    pendingMutations[coord] = list = new List<TileMutation>();
                list.Add(new TileMutation { Coord = coord, PolygonIndex = polygonIndex, Value = value });
            }
        }

        public bool TryGetTileProperty(ChunkCoord coord, int polygonIndex, out TileProperty value)
        {
            lock (tileLock)
            {
                value = default;
                if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return false;
                if (polygonIndex < 0 || polygonIndex >= chunk.Polygons.Count) return false;
                if (chunk.TileProperties != null) value = chunk.TileProperties[polygonIndex];
                return true;
            }
        }

        // Replay any mutations queued against a coord before it was loaded. Called from the
        // promotion paths (EnsureChunk / TryReactivate / AddChunk) right after the chunk lands
        // in `chunks`, so consumers never observe the chunk in an "almost painted" state.
        void ReplayPendingMutations(ChunkCoord coord)
        {
            lock (tileLock)
            {
                if (!pendingMutations.TryGetValue(coord, out List<TileMutation> list)) return;
                pendingMutations.Remove(coord);
                if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return;
                if (chunk.TileProperties == null)
                    chunk.TileProperties = new TileProperty[chunk.Polygons.Count];
                for (int i = 0; i < list.Count; i++)
                {
                    TileMutation m = list[i];
                    if (m.PolygonIndex < 0 || m.PolygonIndex >= chunk.Polygons.Count) continue;
                    chunk.TileProperties[m.PolygonIndex] = m.Value;
                }
                chunk.Version++;
            }
        }

        Vector3 ChunkToWorld(ChunkCoord coord) => coord.WorldCenter(hexRadius, chunkGridSize);

        // Make a chunk active: revive from cache if it was recently evicted, restore from the
        // store if we've persisted it before, otherwise generate it fresh. Returns true only
        // when the chunk was freshly generated (born at flat lattice positions) — revived /
        // restored chunks come back already relaxed, so the caller relaxes only fresh seams.
        public bool EnsureChunk(ChunkCoord coord)
        {
            if (chunks.ContainsKey(coord)) return false;
            if (TryReactivate(coord)) return false;
            chunks[coord] = Generate(coord);
            ReplayPendingMutations(coord);
            return true;
        }

        // Try to make a chunk active without generating. Order: already loaded → cache hit
        // (cheap, in-memory move) → store hit (deserialize). Returns true on any hit. The
        // streamer uses this in classification: anything that fails goes to parallel generation.
        // Cascades on promotion (cache or store): our 6 neighbours' duals may have been built
        // without us while we were away — they'll rebuild on the next BuildDual.
        public bool TryReactivate(ChunkCoord coord)
        {
            if (chunks.ContainsKey(coord)) return true;
            if (cache.TryGetValue(coord, out PrimalChunk cached))
            {
                cache.Remove(coord);
                chunks[coord] = cached;
                InvalidateNeighborDualsOnAdd(coord);
                ReplayPendingMutations(coord);
                return true;
            }
            if (Store != null && Store.Has(coord))
            {
                chunks[coord] = ChunkSnapshot.Restore(coord, Store.Load(coord));
                InvalidateNeighborDualsOnAdd(coord);
                ReplayPendingMutations(coord);
                return true;
            }
            return false;
        }

        // --- Parallel-generation hooks ---
        // Generation is self-contained (own RNG + VertexCollection, pure math, no shared state),
        // so the streamer can produce many fresh chunks across threads and then insert them
        // serially. These expose the decision/insert steps that EnsureChunk normally bundles.

        // True if this coord has a saved snapshot to restore (cheap, serial — touches the store).
        public bool HasStored(ChunkCoord coord) => Store != null && Store.Has(coord);

        // Thread-safe: builds and returns a fresh chunk without touching shared state. Does not
        // insert it — the caller adds it via AddChunk on the owning thread.
        public PrimalChunk GenerateChunk(ChunkCoord coord) => Generate(coord);

        // Insert a generated chunk into the authoritative set (serial). Cascades to neighbours
        // whose duals don't already account for this coord: a neighbour with DualComplete=true
        // was built when all 6 of its neighbours were loaded — including this coord — so its
        // dual is already correct and we skip invalidation. Only DualComplete=false neighbours
        // (built with some coord missing) need rebuilding to incorporate us.
        public void AddChunk(PrimalChunk chunk)
        {
            chunks[chunk.Coord] = chunk;
            InvalidateNeighborDualsOnAdd(chunk.Coord);
            ReplayPendingMutations(chunk.Coord);
        }

        // Generate a chunk's primal grid (interior relaxed, borders pinned at deterministic
        // lattice positions), then bake in the elevation field on every vertex and derive the
        // per-cell terrain classification + chunk-level IsFlat flag from the resulting heights.
        PrimalChunk Generate(ChunkCoord coord)
        {
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.GenerateChunk);

            int chunkSeed = seed ^ (coord.q * 73856093) ^ (coord.r * 19349663);
            var random = new System.Random(chunkSeed);

            List<Polygon> polygons;
            VertexCollection verts;
            using (TerrainProfiler.Measure(TerrainProfiler.Phase.GenPrimal))
                (polygons, verts) = PolygonGridGenerator.GeneratePrimal(
                    chunkGridSize, hexRadius, random, ChunkToWorld(coord), nbIterRelaxation, normalizedRelaxation);

            // 1) Bake elevation. Vertex positions are still at Y=0 at this point — sampling
            //    here makes the primal a true 3D mesh ready for both rendering and gameplay.
            //    Border relax (in RelaxBorders) preserves Y, so this is the only place where Y
            //    is set on a fresh chunk.
            //    IsFlat = "every vertex bottomed out at the bit-exact seabed Y" — anything
            //    above SeabedY (beach band ramp, land profile) means the dual mesh has real
            //    geometry and we cannot short-circuit to the shared flat tile.
            bool isFlat = true;
            float seabedY = Elevation.SeabedY;
            using (TerrainProfiler.Measure(TerrainProfiler.Phase.GenElevation))
            {
                foreach (Vertex v in verts.ToArray())
                {
                    float h = Elevation.Sample(v.Position.x, v.Position.z);
                    v.SetHeight(h);
                    if (h > seabedY) isFlat = false;
                }
            }

            // 2) Classify cells from their vertex heights. Cheap O(verts/cell) per polygon.
            //    "Land" means above the water line (Y=0); the seabed and beach-band slope are
            //    both underwater and classify as ocean — gameplay (NPCs, building placement)
            //    treats those as un-walkable.
            using (TerrainProfiler.Measure(TerrainProfiler.Phase.GenClassify))
            {
                foreach (Polygon p in polygons)
                {
                    bool anyLand = false, anyOcean = false;
                    foreach (Vertex v in p.GetVertices())
                    {
                        if (v.Position.y > 0f) anyLand = true;
                        else anyOcean = true;
                        if (anyLand && anyOcean) break;
                    }
                    p.Terrain = anyLand
                        ? (anyOcean ? CellTerrain.Coastal : CellTerrain.Land)
                        : CellTerrain.Ocean;
                }
            }

            TerrainProfiler.IncChunksGenerated();
            return new PrimalChunk(coord, polygons, verts) { IsFlat = isFlat };
        }

        // Free a chunk from active memory, first saving its (relaxed) state so a later reload
        // restores it exactly. Bypasses the in-memory cache — use SoftUnloadChunk for the
        // normal streaming eviction path.
        public void UnloadChunk(ChunkCoord coord)
        {
            if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return;
            Store?.Save(coord, ChunkSnapshot.Capture(chunk));
            chunks.Remove(coord);
        }

        // Soft eviction: pop the chunk out of the active set into the in-memory cache. The
        // chunk's live state (relaxed vertices, polygon graph) is preserved so the next
        // request can promote it back via TryReactivate — no serialize/deserialize round trip
        // and no store I/O. Cache size is bounded by TrimCache.
        public void SoftUnloadChunk(ChunkCoord coord)
        {
            if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return;
            chunks.Remove(coord);
            cache[coord] = chunk;
            cacheOrder.Enqueue(coord);
        }

        // Flush oldest cached chunks to the store until cache.Count ≤ CacheThreshold. Stale
        // queue entries (chunks that were promoted out before reaching the front) are skipped.
        // Cheap when the cache is under threshold; the streamer calls this once per pass.
        public void TrimCache()
        {
            while (cache.Count > CacheThreshold && cacheOrder.Count > 0)
            {
                ChunkCoord coord = cacheOrder.Dequeue();
                if (!cache.TryGetValue(coord, out PrimalChunk chunk)) continue; // promoted out
                cache.Remove(coord);
                Store?.Save(coord, ChunkSnapshot.Capture(chunk));
            }
        }

        // --- Render-side mirror operations (main thread) ---
        // The streamer owns the authoritative model on a worker thread and publishes
        // snapshots; a second TerrainModel on the main thread mirrors those for the renderer
        // to read. These apply a publish without generating or persisting — the worker is the
        // sole owner of generation and the chunk store.

        // Install a chunk on the main-thread mirror. The chunk arrives with its dual already
        // built by the worker — and the worker also rebuilt any neighbours whose duals our
        // arrival invalidated, so those are coming in this same publish batch. No mirror-
        // side cascade needed.
        public void Install(PrimalChunk chunk) => chunks[chunk.Coord] = chunk;

        // Drop a mirrored chunk without persisting (the worker owns persistence).
        public void Evict(ChunkCoord coord) => chunks.Remove(coord);

        // --- Worker-side dual computation ---
        // After load/relax and before publish, the streamer asks for any chunks whose dual
        // is stale to be (re)built — in parallel. The dual is then carried along with the
        // primal in the deep-copy to consumers.

        // Worker-side check: does this chunk need its dual (re)built? Flat chunks never need
        // one — the surface mounts a shared flat-ocean tile instead of meshing their polygons.
        public bool NeedsDual(ChunkCoord coord)
        {
            return chunks.TryGetValue(coord, out PrimalChunk c)
                && !c.IsFlat
                && (c.Dual == null || c.DualBuiltFromVersion != c.Version);
        }

        // Build this chunk's dual using its own primal + loaded neighbours' primals as
        // cross-chunk context (for completing border cells and computing deterministic
        // ownership). Safe to call in parallel for different coords: reads model.chunks
        // (no concurrent mutation at this point of the pass) and writes only chunk.Dual
        // on the specific coord.
        //
        // Flat chunks short-circuit: we set Dual to a non-null empty sentinel so the
        // version/publish machinery treats their dual as "built" (zero polygons), and the
        // surface uses the IsFlat flag to mount a shared flat tile.
        public void BuildDual(ChunkCoord coord)
        {
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.BuildDual);

            if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return;
            if (chunk.IsFlat)
            {
                chunk.Dual = noPolygons;
                chunk.DualBuiltFromVersion = chunk.Version;
                chunk.DualComplete = true; // flat tiles don't depend on neighbours at all
                TerrainProfiler.IncDualsBuilt();
                return;
            }
            if (chunk.Dual != null && chunk.DualBuiltFromVersion == chunk.Version) return;

            Dictionary<(int, int), List<Polygon>> neighborFaces;
            Dictionary<(int, int), ChunkCoord> minNeighborCoord;

            int loadedNeighbors = 0; // counted during gather; feeds DualComplete at the end
            using (TerrainProfiler.Measure(TerrainProfiler.Phase.BuildDualGather))
            {
                neighborFaces = new Dictionary<(int, int), List<Polygon>>();
                minNeighborCoord = new Dictionary<(int, int), ChunkCoord>();

                foreach (ChunkCoord n in coord.HexesInRange(1))
                {
                    if (n == coord) continue;
                    if (!chunks.TryGetValue(n, out PrimalChunk np)) continue;
                    loadedNeighbors++;

                    // Flat chunks never build their dual (they use a shared flat tile at present
                    // time), so they must not enter the ownership race — otherwise a seam cell
                    // whose lowest-coord chunk is flat would be orphaned (skipped by everyone).
                    // Their primal faces *are* still useful for completing coastal cells on the
                    // non-flat side: a flat neighbour's polygon centers sit at Y=0, which gives
                    // the correct sea-level vertex on a beach dual cell.
                    bool npBuildsDual = !np.IsFlat;

                    foreach (Vertex v in np.Verts.Values)
                    {
                        if (!v.IsEdge) continue;
                        (int, int) key = PolygonGridGenerator.LatticeKey(v.Position, hexRadius);

                        if (!neighborFaces.TryGetValue(key, out List<Polygon> list))
                            neighborFaces[key] = list = new List<Polygon>();
                        foreach (Polygon p in v.Polygons) list.Add(p);

                        if (npBuildsDual
                            && (!minNeighborCoord.TryGetValue(key, out ChunkCoord cur) || n.CompareTo(cur) < 0))
                            minNeighborCoord[key] = n;
                    }
                }
            }

            using (TerrainProfiler.Measure(TerrainProfiler.Phase.BuildDualGenerate))
            {
                var (dualPolygons, _) = PolygonGridGenerator.GenerateDual(
                    chunk, hexRadius, neighborFaces, minNeighborCoord);
                chunk.Dual = dualPolygons;        // GenerateDual now returns List directly — no ToList needed.
                chunk.DualBuiltFromVersion = chunk.Version;
                // "Complete" means every one of the 6 candidate neighbour coords was loaded
                // — i.e. ownership and seam completion were fully resolved. The add-cascade
                // skips invalidating complete neighbours, which makes camera-pan revisits free.
                chunk.DualComplete = loadedNeighbors == 6;
            }
            TerrainProfiler.IncDualsBuilt();
        }

        // Unconditional cascade. Used by RelaxBorders, where the chunk's primal actually moved
        // and any neighbour's dual that referenced its old positions is genuinely stale —
        // DualComplete doesn't save us, the data really has changed.
        void InvalidateNeighborDuals(ChunkCoord coord)
        {
            foreach (ChunkCoord n in coord.HexesInRange(1))
            {
                if (n == coord) continue;
                if (chunks.TryGetValue(n, out PrimalChunk active)) InvalidateDual(active);
                else if (cache.TryGetValue(n, out PrimalChunk cached)) InvalidateDual(cached);
            }
        }

        // Conditional cascade. Used when a chunk is *added* to the loaded set (fresh generate,
        // cache revive, or store restore). A neighbour whose dual is DualComplete was built
        // when all 6 of its candidate neighbour coords were loaded — including this newcomer's
        // coord — so its ownership and seam cells are already correct and we skip it. Only
        // neighbours with DualComplete=false (built with some neighbour missing) need rebuilding
        // to acknowledge the newcomer. This is the optimization that makes camera-pan revisits
        // free: cached chunks come back together, find each other already complete, and no work
        // cascades.
        void InvalidateNeighborDualsOnAdd(ChunkCoord coord)
        {
            foreach (ChunkCoord n in coord.HexesInRange(1))
            {
                if (n == coord) continue;
                if (chunks.TryGetValue(n, out PrimalChunk active))
                {
                    if (!active.DualComplete) InvalidateDual(active);
                }
                else if (cache.TryGetValue(n, out PrimalChunk cached))
                {
                    if (!cached.DualComplete) InvalidateDual(cached);
                }
            }
        }

        static void InvalidateDual(PrimalChunk c)
        {
            if (c.Dual == null) return; // already stale; don't re-bump version
            c.Dual = null;
            c.DualBuiltFromVersion = -1;
            c.DualComplete = false;     // semantic: the (non-existent) dual no longer "knows about" anything
            c.Version++;
        }

        // Joint border relaxation over only the `active` chunks (the freshly-generated frontier
        // plus their loaded neighbours). Restored and settled chunks are deliberately excluded —
        // their seams are already relaxed and welded, so re-relaxing them is pure cost and would
        // bump their versions, forcing needless republish/restore/re-mesh. Including a fresh
        // chunk's loaded neighbours lets the new seam weld to the established (settled) side.
        // Bumps the version of every active chunk whose vertices actually moved.
        public void RelaxBorders(IReadOnlyCollection<ChunkCoord> active)
        {
            if (active == null || active.Count == 0) return;
            using var _ = TerrainProfiler.Measure(TerrainProfiler.Phase.RelaxBorders);
            TerrainProfiler.IncRelaxPasses();

            var coords = new List<ChunkCoord>(active.Count);
            var collections = new List<VertexCollection>(active.Count);
            foreach (ChunkCoord coord in active)
                if (chunks.TryGetValue(coord, out PrimalChunk pc))
                {
                    coords.Add(coord);
                    collections.Add(pc.Verts);
                }

            HashSet<int> moved = PolygonGridGenerator.RelaxBorders(
                collections, hexRadius, nbIterBorderRelaxation, borderRelaxInteriorRings);

            foreach (int i in moved)
            {
                PrimalChunk c = chunks[coords[i]];
                c.Version++;
                // Own dual is stale (vertices moved); kill it directly without re-bumping the
                // version we just incremented.
                c.Dual = null;
                c.DualBuiltFromVersion = -1;
                // Cascade: neighbours' seam cells used our pre-move quad positions.
                InvalidateNeighborDuals(c.Coord);
            }
        }
    }
}
