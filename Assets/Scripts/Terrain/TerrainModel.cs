using System.Collections.Generic;
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

        // Soft-eviction tier: chunks that no consumer wants right now but are kept in memory
        // so the next pass can promote them back with no I/O. Promotion happens via
        // TryReactivate; overflow flushes to Store via TrimCache. cacheOrder is FIFO and may
        // contain stale entries (promoted-out coords), which are skipped during trimming.
        readonly Dictionary<ChunkCoord, PrimalChunk> cache = new();
        readonly Queue<ChunkCoord> cacheOrder = new();

        public IEnumerable<ChunkCoord> LoadedCoords => chunks.Keys;
        public int LoadedCount => chunks.Count;
        public int CachedCount => cache.Count;
        public bool IsLoaded(ChunkCoord coord) => chunks.ContainsKey(coord);
        public bool IsCached(ChunkCoord coord) => cache.ContainsKey(coord);
        public bool TryGet(ChunkCoord coord, out PrimalChunk chunk) => chunks.TryGetValue(coord, out chunk);

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
            return true;
        }

        // Try to make a chunk active without generating. Order: already loaded → cache hit
        // (cheap, in-memory move) → store hit (deserialize). Returns true on any hit. The
        // streamer uses this in classification: anything that fails goes to parallel generation.
        public bool TryReactivate(ChunkCoord coord)
        {
            if (chunks.ContainsKey(coord)) return true;
            if (cache.TryGetValue(coord, out PrimalChunk cached))
            {
                cache.Remove(coord);
                chunks[coord] = cached;
                return true;
            }
            if (Store != null && Store.Has(coord))
            {
                chunks[coord] = ChunkSnapshot.Restore(coord, Store.Load(coord));
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

        // Insert a generated chunk into the authoritative set (serial).
        public void AddChunk(PrimalChunk chunk) => chunks[chunk.Coord] = chunk;

        // Generate a chunk's primal grid (interior relaxed, borders pinned at deterministic
        // lattice positions).
        PrimalChunk Generate(ChunkCoord coord)
        {
            int chunkSeed = seed ^ (coord.q * 73856093) ^ (coord.r * 19349663);
            var random = new System.Random(chunkSeed);

            var (polygons, verts) = PolygonGridGenerator.GeneratePrimal(
                chunkGridSize, hexRadius, random, ChunkToWorld(coord), nbIterRelaxation, normalizedRelaxation);

            return new PrimalChunk(coord, polygons, verts);
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

        // Install a chunk restored from a worker snapshot, replacing any existing one, and
        // drop the cached Dual on each of this chunk's 6 loaded neighbours. We have to
        // cascade in two cases:
        //   • New chunk: neighbours that built their dual before this chunk arrived had
        //     LookupPrimal return null for our coord, so their X-Y seam cells were dropped
        //     (<3 faces). They need to rebuild now that we exist.
        //   • New version: neighbours' seam cells were built against the previous version's
        //     quad positions and are now stale.
        // Both cases are handled by unconditionally invalidating neighbours. They rebuild
        // lazily on the next ChunkSurface.Apply().
        public void Install(PrimalChunk chunk)
        {
            foreach (ChunkCoord n in chunk.Coord.HexesInRange(1))
                if (n != chunk.Coord && chunks.TryGetValue(n, out PrimalChunk neighbor))
                {
                    neighbor.Dual = null;
                    neighbor.DualBuiltFromVersion = -1;
                }
            chunks[chunk.Coord] = chunk;
        }

        // Drop a mirrored chunk without persisting (the worker owns persistence).
        public void Evict(ChunkCoord coord) => chunks.Remove(coord);

        // Joint border relaxation over only the `active` chunks (the freshly-generated frontier
        // plus their loaded neighbours). Restored and settled chunks are deliberately excluded —
        // their seams are already relaxed and welded, so re-relaxing them is pure cost and would
        // bump their versions, forcing needless republish/restore/re-mesh. Including a fresh
        // chunk's loaded neighbours lets the new seam weld to the established (settled) side.
        // Bumps the version of every active chunk whose vertices actually moved.
        public void RelaxBorders(IReadOnlyCollection<ChunkCoord> active)
        {
            if (active == null || active.Count == 0) return;

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
                chunks[coords[i]].Version++;
        }
    }
}
