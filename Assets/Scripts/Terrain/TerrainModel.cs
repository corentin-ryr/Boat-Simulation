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
                    // Cross-chunk cascade. The road mesh on a neighbour chunk reads our
                    // TileProperties to decide whether to extend an arm across the seam, so
                    // if the painted cell touches a seam we must bump the neighbour chunks'
                    // Versions to retrigger BuildRoadMesh on that side. Without this, the
                    // first paint on each side of a seam wouldn't visually connect until
                    // some unrelated change bumped the neighbour.
                    BumpNeighborChunksForCellEdges(chunk, coord, polygonIndex);
                    return;
                }

                if (!pendingMutations.TryGetValue(coord, out List<TileMutation> list))
                    pendingMutations[coord] = list = new List<TileMutation>();
                list.Add(new TileMutation { Coord = coord, PolygonIndex = polygonIndex, Value = value });
            }
        }

        // For each cross-chunk edge of the painted cell, bump that neighbour chunk's
        // Version. Precise path: walk the chunk's per-edge neighbour table built by
        // BuildPrimalNeighborTable. Fallback (table null — e.g. a fresh chunk being
        // painted before its first dual build): if any of the cell's primal vertices is on
        // the chunk boundary, conservatively bump all 6 neighbour chunks. The fallback over-
        // bumps but is correct; the precise path runs in steady state.
        void BumpNeighborChunksForCellEdges(PrimalChunk chunk, ChunkCoord coord, int polygonIndex)
        {
            int[] offsets = chunk.PrimalEdgeOffsets;
            int[] dqArr   = chunk.PrimalNeighborChunkDQ;
            int[] drArr   = chunk.PrimalNeighborChunkDR;
            int[] idxArr  = chunk.PrimalNeighborPolyIdx;

            if (offsets != null && dqArr != null && drArr != null && idxArr != null
                && polygonIndex + 1 < offsets.Length)
            {
                int eStart = offsets[polygonIndex];
                int eEnd = offsets[polygonIndex + 1];
                // Use a small bitset to coalesce repeated neighbour bumps within the cell.
                // 6 neighbour chunks max — but tracking by (dq,dr) means we just dedupe
                // through an inline list scan, no hashing needed.
                int bumpCount = 0;
                System.Span<int> dqBuf = stackalloc int[6];
                System.Span<int> drBuf = stackalloc int[6];
                for (int j = eStart; j < eEnd; j++)
                {
                    int dq = dqArr[j], dr = drArr[j];
                    if (dq == 0 && dr == 0) continue;       // own-chunk edge
                    if (idxArr[j] < 0) continue;            // unresolved seam
                    bool dup = false;
                    for (int k = 0; k < bumpCount; k++)
                        if (dqBuf[k] == dq && drBuf[k] == dr) { dup = true; break; }
                    if (dup) continue;
                    if (bumpCount < dqBuf.Length) { dqBuf[bumpCount] = dq; drBuf[bumpCount] = dr; bumpCount++; }
                    ChunkCoord nc = new ChunkCoord(coord.q + dq, coord.r + dr);
                    if (chunks.TryGetValue(nc, out PrimalChunk nbChunk))
                        nbChunk.Version++;
                }
                return;
            }

            // Fallback path: table not built yet. Bump all 6 neighbour chunks if the cell
            // touches the chunk boundary (any vertex IsEdge). Cheap to check and correct
            // in all cases — just less precise.
            Polygon poly = chunk.Polygons != null && polygonIndex < chunk.Polygons.Count
                ? chunk.Polygons[polygonIndex] : null;
            if (poly == null) return;
            bool onSeam = false;
            foreach (Vertex v in poly.GetVertices()) if (v.IsEdge) { onSeam = true; break; }
            if (!onSeam) return;
            foreach (ChunkCoord nc in coord.HexesInRange(1))
            {
                if (nc == coord) continue;
                if (chunks.TryGetValue(nc, out PrimalChunk nbChunk)) nbChunk.Version++;
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
                    // Cascade per painted cell so neighbour chunks rebuild their road mesh
                    // and pick up cross-seam connections introduced by the replayed paints.
                    BumpNeighborChunksForCellEdges(chunk, coord, m.PolygonIndex);
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
        //
        // Also (re)builds the chunk's primal neighbour table while we hold full neighbour
        // context: store-restored chunks arrived without one, and cache-revived chunks had
        // theirs frozen at soft-unload time so it may reference neighbours that have come or
        // gone in the interim.
        public bool TryReactivate(ChunkCoord coord)
        {
            if (chunks.ContainsKey(coord)) return true;
            if (cache.TryGetValue(coord, out PrimalChunk cached))
            {
                cache.Remove(coord);
                chunks[coord] = cached;
                InvalidateNeighborDualsOnAdd(coord);
                BuildPrimalNeighborTable(coord);
                ReplayPendingMutations(coord);
                return true;
            }
            if (Store != null && Store.Has(coord))
            {
                chunks[coord] = ChunkSnapshot.Restore(coord, Store.Load(coord));
                InvalidateNeighborDualsOnAdd(coord);
                BuildPrimalNeighborTable(coord);
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
            // Build the chunk's primal neighbour table now that we hold the model lock and
            // can see whatever neighbour chunks were loaded earlier in this same pass.
            // Cross-chunk slots whose neighbour isn't loaded yet stay at -1; the dual cascade
            // will rebuild this when a missing neighbour shows up (BuildDual end).
            BuildPrimalNeighborTable(chunk.Coord);
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

            // Cache primal cell centroids so the render-side painter can do nearest-centroid
            // lookups without walking the primal graph (which DeepCopy drops). Cheap O(verts)
            // per polygon; shipped through DeepCopy onto the consumer mirror.
            Vector3[] centroids = new Vector3[polygons.Count];
            for (int i = 0; i < polygons.Count; i++)
                centroids[i] = polygons[i].GetCenter();

            TerrainProfiler.IncChunksGenerated();
            return new PrimalChunk(coord, polygons, verts)
            {
                IsFlat = isFlat,
                PrimalCellCentroids = centroids,
            };
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
                // Reverse lookup: "what's the index of this Polygon in chunk.Polygons?"
                // Fed to GenerateDual so each dual corner can be tagged with the owning
                // primal cell index (used by the render-side vertex-colour tint). Polygons
                // from neighbouring chunks aren't in this dict and naturally fall to -1
                // (= neutral palette entry at render time).
                var primalIndex = new Dictionary<Polygon, int>(chunk.Polygons.Count);
                for (int i = 0; i < chunk.Polygons.Count; i++) primalIndex[chunk.Polygons[i]] = i;

                var (dualPolygons, _, dualVertPrimalIdx) = PolygonGridGenerator.GenerateDual(
                    chunk, hexRadius, neighborFaces, minNeighborCoord, primalIndex);
                chunk.Dual = dualPolygons;        // GenerateDual now returns List directly — no ToList needed.
                chunk.DualVertexPrimalIndex = dualVertPrimalIdx;
                chunk.DualBuiltFromVersion = chunk.Version;
                // "Complete" means every one of the 6 candidate neighbour coords was loaded
                // — i.e. ownership and seam completion were fully resolved. The add-cascade
                // skips invalidating complete neighbours, which makes camera-pan revisits free.
                chunk.DualComplete = loadedNeighbors == 6;
            }

            // Same context (loaded neighbour set) is exactly what BuildPrimalNeighborTable
            // needs, so refresh the chunk's per-edge neighbour table now. If a neighbour
            // shows up later, the dual cascade reinvalidates this chunk and brings us back
            // here — keeping the table in sync with the dual.
            BuildPrimalNeighborTable(coord);

            TerrainProfiler.IncDualsBuilt();
        }

        // Build the per-cell, per-edge neighbour table for the chunk. Each edge slot stores
        // its two endpoint positions plus a reference to the cell across that edge: (0,0,i)
        // for an in-chunk neighbour, (dq,dr,i) for a cell in a loaded neighbour chunk, or
        // (0,0,-1) when the neighbour cell can't be resolved (chunk seam where no neighbour
        // chunk is loaded). Consumed by ApplyTileMutation (cross-chunk Version cascade on
        // paint) and by ChunkSurface.BuildRoadMesh (road arms connect through this).
        //
        // Worker-thread only; reads the model's `chunks` dict to resolve cross-chunk
        // neighbours, so safe to call from anywhere the BuildDual pass runs (RunPass holds
        // the world steady around these calls). Flat chunks skip — they don't paint roads.
        public void BuildPrimalNeighborTable(ChunkCoord coord)
        {
            if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return;
            if (chunk.IsFlat) return;
            List<Polygon> polygons = chunk.Polygons;
            if (polygons == null) return;

            int cellCount = polygons.Count;
            int totalEdges = 0;
            for (int i = 0; i < cellCount; i++) totalEdges += polygons[i].GetVertices().Length;

            int[]     offsets = new int[cellCount + 1];
            Vector3[] vAarr  = new Vector3[totalEdges];
            Vector3[] vBarr  = new Vector3[totalEdges];
            int[]     dqArr  = new int[totalEdges];
            int[]     drArr  = new int[totalEdges];
            int[]     idxArr = new int[totalEdges];

            // Own-chunk reverse map: Polygon → index. Built once, used to convert
            // p.NeighborAtEdge(j) into an integer index when the neighbour is in-chunk.
            var ownIndex = new Dictionary<Polygon, int>(cellCount);
            for (int i = 0; i < cellCount; i++) ownIndex[polygons[i]] = i;

            // Per-neighbour-chunk reverse edge map: (orderedLatticeKeyA, orderedLatticeKeyB)
            // → polygon index in that chunk. Built once per neighbour, then probed per
            // cross-chunk edge of the painted chunk. Cheap: ~50 cells × ~4 edges per
            // neighbour, 6 neighbours, all dictionary work.
            var neighborEdgeMaps = new Dictionary<ChunkCoord,
                Dictionary<((int, int), (int, int)), int>>();
            foreach (ChunkCoord nc in coord.HexesInRange(1))
            {
                if (nc == coord) continue;
                if (!chunks.TryGetValue(nc, out PrimalChunk np)) continue;
                if (np.IsFlat) continue;

                var edgeMap = new Dictionary<((int, int), (int, int)), int>();
                List<Polygon> npPolys = np.Polygons;
                for (int pi = 0; pi < npPolys.Count; pi++)
                {
                    Vertex[] vs = npPolys[pi].GetVertices();
                    int n = vs.Length;
                    for (int j = 0; j < n; j++)
                    {
                        Vertex va = vs[j];
                        Vertex vb = vs[(j + 1) % n];
                        // Only edges with BOTH endpoints on the chunk boundary can lie on
                        // a cross-chunk seam — interior edges can't match anything in
                        // another chunk and would just clutter the map.
                        if (!va.IsEdge || !vb.IsEdge) continue;
                        (int, int) ka = PolygonGridGenerator.LatticeKey(va.Position, hexRadius);
                        (int, int) kb = PolygonGridGenerator.LatticeKey(vb.Position, hexRadius);
                        edgeMap[OrderedEdgeKey(ka, kb)] = pi;
                    }
                }
                neighborEdgeMaps[nc] = edgeMap;
            }

            int cursor = 0;
            for (int i = 0; i < cellCount; i++)
            {
                offsets[i] = cursor;
                Polygon p = polygons[i];
                Vertex[] verts = p.GetVertices();
                int n = verts.Length;

                for (int j = 0; j < n; j++)
                {
                    Vertex va = verts[j];
                    Vertex vb = verts[(j + 1) % n];
                    vAarr[cursor] = va.Position;
                    vBarr[cursor] = vb.Position;

                    // In-chunk neighbour across edge (va, vb): the unique polygon that
                    // shares BOTH endpoints with `p`. Walk vertex.Polygons (always
                    // populated by the Polygon ctor, on both fresh-generated and
                    // snapshot-restored chunks) instead of p.GetNeighbors() — the latter
                    // is set by ComputeNeighbors which Restore doesn't call.
                    Polygon nb = null;
                    foreach (Polygon cand in va.Polygons)
                    {
                        if (cand == p) continue;
                        bool sharesVb = false;
                        foreach (Polygon other in vb.Polygons)
                            if (other == cand) { sharesVb = true; break; }
                        if (sharesVb) { nb = cand; break; }
                    }
                    if (nb != null && ownIndex.TryGetValue(nb, out int nbIdx))
                    {
                        dqArr[cursor]  = 0;
                        drArr[cursor]  = 0;
                        idxArr[cursor] = nbIdx;
                    }
                    else
                    {
                        // Cross-chunk: resolve via lattice keys on the two endpoints. Look
                        // in every loaded neighbour-chunk edge map for an entry with the
                        // exact same undirected key.
                        int foundIdx = -1;
                        ChunkCoord foundCoord = default;
                        if (va.IsEdge && vb.IsEdge)
                        {
                            (int, int) ka = PolygonGridGenerator.LatticeKey(va.Position, hexRadius);
                            (int, int) kb = PolygonGridGenerator.LatticeKey(vb.Position, hexRadius);
                            var key = OrderedEdgeKey(ka, kb);
                            foreach (var kv in neighborEdgeMaps)
                            {
                                if (kv.Value.TryGetValue(key, out int hit))
                                {
                                    foundIdx = hit;
                                    foundCoord = kv.Key;
                                    break;
                                }
                            }
                        }

                        if (foundIdx >= 0)
                        {
                            dqArr[cursor]  = foundCoord.q - coord.q;
                            drArr[cursor]  = foundCoord.r - coord.r;
                            idxArr[cursor] = foundIdx;
                        }
                        else
                        {
                            // No neighbour loaded right now (or interior boundary).
                            // BuildPrimalNeighborTable will be re-run when something
                            // changes (dual cascade, neighbour add) and may resolve.
                            dqArr[cursor]  = 0;
                            drArr[cursor]  = 0;
                            idxArr[cursor] = -1;
                        }
                    }
                    cursor++;
                }
            }
            offsets[cellCount] = cursor;

            chunk.PrimalEdgeOffsets      = offsets;
            chunk.PrimalEdgeVertexA      = vAarr;
            chunk.PrimalEdgeVertexB      = vBarr;
            chunk.PrimalNeighborChunkDQ  = dqArr;
            chunk.PrimalNeighborChunkDR  = drArr;
            chunk.PrimalNeighborPolyIdx  = idxArr;
        }

        // Order-independent key for an undirected edge whose endpoints are lattice keys —
        // mirrors PolygonGridGenerator.EdgeKey but on rounded integer coords, so a chunk's
        // seam edge produces the SAME key as the matching edge from the chunk across the
        // seam regardless of which winding each side traversed it in.
        static ((int, int), (int, int)) OrderedEdgeKey((int, int) a, (int, int) b)
        {
            bool aFirst = a.Item1 < b.Item1 || (a.Item1 == b.Item1 && a.Item2 < b.Item2);
            return aFirst ? (a, b) : (b, a);
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
