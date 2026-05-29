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
        // Cascades on promotion (cache or store): our 6 neighbours' duals may have been built
        // without us while we were away — they'll rebuild on the next BuildDual.
        public bool TryReactivate(ChunkCoord coord)
        {
            if (chunks.ContainsKey(coord)) return true;
            if (cache.TryGetValue(coord, out PrimalChunk cached))
            {
                cache.Remove(coord);
                chunks[coord] = cached;
                InvalidateNeighborDuals(coord);
                return true;
            }
            if (Store != null && Store.Has(coord))
            {
                chunks[coord] = ChunkSnapshot.Restore(coord, Store.Load(coord));
                InvalidateNeighborDuals(coord);
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

        // Insert a generated chunk into the authoritative set (serial). Cascades: this chunk's
        // 6 loaded-or-cached neighbours had their cached duals (if any) built without us — now
        // stale. Invalidate them so the next BuildDual rebuilds with our data.
        public void AddChunk(PrimalChunk chunk)
        {
            chunks[chunk.Coord] = chunk;
            InvalidateNeighborDuals(chunk.Coord);
        }

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

        // Worker-side check: does this chunk need its dual (re)built?
        public bool NeedsDual(ChunkCoord coord)
        {
            return chunks.TryGetValue(coord, out PrimalChunk c)
                && (c.Dual == null || c.DualBuiltFromVersion != c.Version);
        }

        // Build this chunk's dual using its own primal + loaded neighbours' primals as
        // cross-chunk context (for completing border cells and computing deterministic
        // ownership). Safe to call in parallel for different coords: reads model.chunks
        // (no concurrent mutation at this point of the pass) and writes only chunk.Dual
        // on the specific coord.
        public void BuildDual(ChunkCoord coord)
        {
            if (!chunks.TryGetValue(coord, out PrimalChunk chunk)) return;
            if (chunk.Dual != null && chunk.DualBuiltFromVersion == chunk.Version) return;

            Dictionary<(int, int), List<Polygon>> neighborFaces = new();
            Dictionary<(int, int), ChunkCoord> minNeighborCoord = new();

            foreach (ChunkCoord n in coord.HexesInRange(1))
            {
                if (n == coord) continue;
                if (!chunks.TryGetValue(n, out PrimalChunk np)) continue;

                foreach (Vertex v in np.Verts.ToArray())
                {
                    if (!v.IsEdge) continue;
                    (int, int) key = PolygonGridGenerator.LatticeKey(v.Position, hexRadius);

                    if (!neighborFaces.TryGetValue(key, out List<Polygon> list))
                        neighborFaces[key] = list = new List<Polygon>();
                    list.AddRange(v.Polygons);

                    if (!minNeighborCoord.TryGetValue(key, out ChunkCoord cur) || n.CompareTo(cur) < 0)
                        minNeighborCoord[key] = n;
                }
            }

            var (dualPolygons, _) = PolygonGridGenerator.GenerateDual(
                chunk, hexRadius, neighborFaces, minNeighborCoord);
            chunk.Dual = dualPolygons.ToList();
            chunk.DualBuiltFromVersion = chunk.Version;
        }

        // Invalidate the cached dual on every loaded-or-cached neighbour of `coord`. The
        // version bump in InvalidateDual makes the publish loop notice — the chunk's primal
        // didn't change but its visible state did. Idempotent: re-invalidating an already-
        // invalidated chunk is a no-op (Dual already null, version not re-bumped).
        void InvalidateNeighborDuals(ChunkCoord coord)
        {
            foreach (ChunkCoord n in coord.HexesInRange(1))
            {
                if (n == coord) continue;
                if (chunks.TryGetValue(n, out PrimalChunk active)) InvalidateDual(active);
                else if (cache.TryGetValue(n, out PrimalChunk cached)) InvalidateDual(cached);
            }
        }

        static void InvalidateDual(PrimalChunk c)
        {
            if (c.Dual == null) return; // already stale; don't re-bump version
            c.Dual = null;
            c.DualBuiltFromVersion = -1;
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
