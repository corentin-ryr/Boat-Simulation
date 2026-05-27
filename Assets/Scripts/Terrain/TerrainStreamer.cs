using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace TerrainGrid
{
    // Handle a consumer holds to declare its desired chunks and drain ready ones. Opaque —
    // obtained from TerrainStreamer.RegisterConsumer.
    public sealed class ChunkConsumer
    {
        internal readonly int Id;
        internal ChunkConsumer(int id) { Id = id; }
    }

    // A completed publish for one consumer. Immutable, thread-safe payload: deep-copied chunks
    // and the consumer's effective set (so it can evict its mirror down to what it still wants).
    public class StreamResult
    {
        public HashSet<ChunkCoord> Loaded;   // the consumer's effective set this pass
        public List<PrimalChunk> Changed;    // independent copies of chunks whose version moved
        public Exception Error;              // worker exception, surfaced to the consumer

        public static StreamResult FromError(Exception e) => new StreamResult { Error = e };
    }

    // Drives the primal data layer on one dedicated background thread, serving any number of
    // consumers. Each consumer declares a desired set of chunk coords and a service level:
    //   • render-ready — chunks are relaxed and surrounded by a 1-ring halo (for seam welding
    //     and the renderer's cross-chunk neighbour lookups), then published as deep copies.
    //   • primal-only  — chunks are loaded as generated (borders pinned), no halo, no relaxation;
    //     isolated/non-contiguous requests are valid. For gameplay/occupancy data.
    // The worker loads the UNION of all consumers' (effective) sets, so a chunk stays alive while
    // any consumer wants it, and is unloaded (persisted to the store) only when none do. Results
    // are routed per consumer. The worker is the sole mutator of the authoritative model.
    public class TerrainStreamer
    {
        class ConsumerState
        {
            public HashSet<ChunkCoord> Desired = new();
            public bool RenderReady;
            public readonly ConcurrentQueue<StreamResult> Results = new();
            // Version last published to THIS consumer per coord — so a chunk is re-copied only
            // when newly wanted or actually changed, independently for each consumer.
            public readonly Dictionary<ChunkCoord, int> PublishedVersion = new();
        }

        readonly TerrainModel model;

        readonly object gate = new object();
        readonly AutoResetEvent signal = new AutoResetEvent(false);
        readonly Dictionary<int, ConsumerState> consumers = new();
        int nextId;
        bool dirty;
        bool stop;
        Thread thread;

        public TerrainStreamer(TerrainModel model) { this.model = model; }

        // --- Consumer API (main thread) ---

        public ChunkConsumer RegisterConsumer(bool renderReady)
        {
            lock (gate)
            {
                int id = nextId++;
                consumers[id] = new ConsumerState { RenderReady = renderReady };
                return new ChunkConsumer(id);
            }
        }

        public void UnregisterConsumer(ChunkConsumer consumer)
        {
            lock (gate) { consumers.Remove(consumer.Id); dirty = true; }
            signal.Set();
        }

        // Declare the consumer's full current desired set (declarative — newest call wins).
        public void SetDesired(ChunkConsumer consumer, IReadOnlyCollection<ChunkCoord> coords)
        {
            lock (gate)
            {
                if (consumers.TryGetValue(consumer.Id, out ConsumerState s))
                {
                    s.Desired = new HashSet<ChunkCoord>(coords);
                    dirty = true;
                }
            }
            signal.Set();
        }

        // Pop a finished publish for this consumer, if any.
        public bool TryDrainResult(ChunkConsumer consumer, out StreamResult result)
        {
            result = null;
            ConsumerState s;
            lock (gate) { consumers.TryGetValue(consumer.Id, out s); }
            return s != null && s.Results.TryDequeue(out result);
        }

        // Run one pass synchronously before the worker starts, so the opening view isn't empty.
        public void Prime() => RunPass();

        public void Start()
        {
            thread = new Thread(WorkerLoop) { Name = "TerrainStreamer", IsBackground = true };
            thread.Start();
        }

        public void Stop()
        {
            lock (gate) { stop = true; }
            signal.Set();
            thread?.Join();
        }

        void WorkerLoop()
        {
            while (true)
            {
                bool has;
                lock (gate) { if (stop) return; has = dirty; }
                if (!has) { signal.WaitOne(); continue; }
                try { RunPass(); }
                catch (Exception e) { EnqueueErrorToAll(e); }
            }
        }

        // True if any consumer changed its request since we latched this pass.
        bool NewerPending() { lock (gate) return dirty; }

        void EnqueueErrorToAll(Exception e)
        {
            lock (gate)
                foreach (ConsumerState s in consumers.Values)
                    s.Results.Enqueue(StreamResult.FromError(e));
        }

        // A set expanded by one hex ring (HexesInRange(1) includes the centre), i.e. set ∪ halo.
        static HashSet<ChunkCoord> ExpandRing(IEnumerable<ChunkCoord> set)
        {
            var result = new HashSet<ChunkCoord>();
            foreach (ChunkCoord c in set)
                foreach (ChunkCoord n in c.HexesInRange(1))
                    result.Add(n);
            return result;
        }

        void RunPass()
        {
            // Snapshot consumer requests under the lock; do all heavy work outside it.
            List<ConsumerState> snap;
            List<HashSet<ChunkCoord>> desiredSnap;
            List<bool> renderSnap;
            lock (gate)
            {
                dirty = false;
                snap = consumers.Values.ToList();
                desiredSnap = snap.Select(s => new HashSet<ChunkCoord>(s.Desired)).ToList();
                renderSnap = snap.Select(s => s.RenderReady).ToList();
            }

            // Loaded set: render consumers' sets expanded by a 1-ring halo (seam welding +
            // neighbour lookups), plus data consumers' raw sets.
            var renderUnion = new HashSet<ChunkCoord>();
            var dataUnion = new HashSet<ChunkCoord>();
            for (int i = 0; i < snap.Count; i++)
                (renderSnap[i] ? renderUnion : dataUnion).UnionWith(desiredSnap[i]);
            HashSet<ChunkCoord> renderLoaded = ExpandRing(renderUnion);
            var loaded = new HashSet<ChunkCoord>(renderLoaded);
            loaded.UnionWith(dataUnion);

            // Evict whatever no consumer needs (persists to the store).
            foreach (ChunkCoord coord in model.LoadedCoords.Where(c => !loaded.Contains(c)).ToList())
                model.UnloadChunk(coord);

            // Classify: restore stored chunks serially; queue the rest for parallel generation.
            var toGenerate = new List<ChunkCoord>();
            foreach (ChunkCoord coord in loaded)
            {
                if (model.IsLoaded(coord)) continue;
                if (model.HasStored(coord)) model.EnsureChunk(coord);
                else toGenerate.Add(coord);
            }

            if (NewerPending()) return; // stale: abandon before the expensive batch

            var fresh = new List<ChunkCoord>();
            if (toGenerate.Count > 0)
            {
                var built = new PrimalChunk[toGenerate.Count];
                Parallel.For(0, toGenerate.Count, k => built[k] = model.GenerateChunk(toGenerate[k]));
                foreach (PrimalChunk pc in built) { model.AddChunk(pc); fresh.Add(pc.Coord); }
            }

            if (NewerPending()) return; // stale: abandon (generated chunks stay loaded for reuse)

            // Relax only fresh chunks that need welding (those in the render-loaded region) plus
            // their loaded neighbours. Data-only chunks keep their pinned borders.
            var active = new HashSet<ChunkCoord>();
            foreach (ChunkCoord f in fresh)
            {
                if (!renderLoaded.Contains(f)) continue;
                active.Add(f);
                foreach (ChunkCoord n in f.HexesInRange(1))
                    if (n != f && model.IsLoaded(n)) active.Add(n);
            }
            model.RelaxBorders(active);

            // Publish per consumer: deep-copy chunks in its effective set whose version changed.
            for (int i = 0; i < snap.Count; i++)
            {
                ConsumerState s = snap[i];
                HashSet<ChunkCoord> effective = renderSnap[i] ? ExpandRing(desiredSnap[i]) : desiredSnap[i];

                var changed = new List<PrimalChunk>();
                foreach (ChunkCoord coord in effective)
                {
                    if (!model.TryGet(coord, out PrimalChunk pc)) continue;
                    if (!s.PublishedVersion.TryGetValue(coord, out int v) || v != pc.Version)
                    {
                        changed.Add(pc.DeepCopy());
                        s.PublishedVersion[coord] = pc.Version;
                    }
                }
                foreach (ChunkCoord stale in s.PublishedVersion.Keys.Where(k => !effective.Contains(k)).ToList())
                    s.PublishedVersion.Remove(stale);

                s.Results.Enqueue(new StreamResult { Loaded = effective, Changed = changed });
            }
        }
    }
}
