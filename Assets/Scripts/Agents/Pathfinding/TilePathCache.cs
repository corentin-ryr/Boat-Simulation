using System.Collections.Generic;
using TerrainGrid;

namespace Agents.Pathfinding
{
    // Bounded cache for TilePathfinder results. Keyed by the (from, to)
    // tile pair. Validity is version-stamped: an entry records the source
    // and target chunks' PrimalChunk.Version at cache time, and the get
    // path re-reads the current versions before returning. Mismatch ⇒
    // the entry is dropped and the caller re-plans.
    //
    // No event subscription needed — version-checking on read is
    // self-correcting. For paths crossing intermediate chunks that change,
    // the agent's IsValid() during MoveTo execution catches the stale
    // step before it walks into an impassable tile (Phase F adds that
    // recheck).
    //
    // FIFO eviction. Capacity 256 is loose at v1 scale; bump or swap to
    // LRU if profile shows cache thrash.
    public sealed class TilePathCache
    {
        readonly Dictionary<Key, Entry> byKey = new Dictionary<Key, Entry>();
        readonly Queue<Key> insertionOrder = new Queue<Key>();
        readonly int capacity;
        readonly ChunkManager mgr;

        public TilePathCache(ChunkManager mgr, int capacity = 256)
        {
            this.mgr = mgr;
            this.capacity = capacity;
        }

        public TilePath GetOrFind(ChunkCoord fromCoord, int fromCell,
                                  ChunkCoord toCoord, int toCell)
        {
            Key key = new Key(fromCoord, fromCell, toCoord, toCell);

            int fromV = ChunkVersionOrZero(fromCoord);
            int toV = ChunkVersionOrZero(toCoord);

            if (byKey.TryGetValue(key, out Entry hit))
            {
                if (hit.FromVersion == fromV && hit.ToVersion == toV)
                    return hit.Path;
                byKey.Remove(key);
                // Stale entries leak into insertionOrder; they're removed
                // when reached by FIFO eviction. No mid-queue removal —
                // keeps the queue a constant-time append.
            }

            TilePath path = TilePathfinder.FindPath(fromCoord, fromCell, toCoord, toCell, mgr);
            if (path == null) return null;

            Put(key, new Entry(path, fromV, toV));
            return path;
        }

        void Put(Key key, Entry entry)
        {
            byKey[key] = entry;
            insertionOrder.Enqueue(key);
            // Evict the oldest insertion that's still live, until we're at
            // or below capacity. Tombstones (stale keys still in the queue
            // but already removed from byKey by version invalidation) are
            // skipped cheaply.
            while (byKey.Count > capacity && insertionOrder.Count > 0)
            {
                Key oldest = insertionOrder.Dequeue();
                byKey.Remove(oldest);
            }
        }

        public void Clear()
        {
            byKey.Clear();
            insertionOrder.Clear();
        }

        public int Count => byKey.Count;

        int ChunkVersionOrZero(ChunkCoord c)
        {
            if (mgr == null) return 0;
            return mgr.TryGetPublishedChunk(c, out PrimalChunk pc) ? pc.Version : -1;
        }

        readonly struct Key : System.IEquatable<Key>
        {
            public readonly ChunkCoord FromCoord;
            public readonly int FromCell;
            public readonly ChunkCoord ToCoord;
            public readonly int ToCell;
            public Key(ChunkCoord fc, int fi, ChunkCoord tc, int ti)
            {
                FromCoord = fc; FromCell = fi; ToCoord = tc; ToCell = ti;
            }
            public bool Equals(Key o)
                => FromCoord.Equals(o.FromCoord) && FromCell == o.FromCell
                && ToCoord.Equals(o.ToCoord) && ToCell == o.ToCell;
            public override bool Equals(object o) => o is Key k && Equals(k);
            public override int GetHashCode()
            {
                unchecked
                {
                    int h = 17;
                    h = h * 31 + FromCoord.GetHashCode();
                    h = h * 31 + FromCell;
                    h = h * 31 + ToCoord.GetHashCode();
                    h = h * 31 + ToCell;
                    return h;
                }
            }
        }

        readonly struct Entry
        {
            public readonly TilePath Path;
            public readonly int FromVersion;
            public readonly int ToVersion;
            public Entry(TilePath p, int fv, int tv) { Path = p; FromVersion = fv; ToVersion = tv; }
        }
    }
}
