using System;
using System.Diagnostics;
using System.Threading;
using Unity.Profiling;

namespace TerrainGrid
{
    // Central profiling for the terrain pipeline. Three complementary views of the same data:
    //
    //   1. ProfilerMarker — shows up in Unity's Profiler window timeline on the appropriate
    //      thread. Use this for per-call inspection (which BuildDual call exactly was slow).
    //
    //   2. Cumulative ms per Phase — the overlay subtracts last-frame from current to read
    //      "how many ms did this phase spend during the frame that just ended", then tracks
    //      the peak over a rolling window. Tells you *which* phase eats the budget on the
    //      worst frame — exactly the question camera-pan slowdowns raise.
    //
    //   3. Throughput counters — chunks generated/built/copied, mesh rebuilds, collider
    //      assignments. Tells you "what's churning, and how hard?"
    //
    // All measurements use Stopwatch ticks (high-resolution, monotonic) and Interlocked
    // updates, so worker-thread phases (BuildDual, RelaxBorders, GenerateChunk) and main-
    // thread phases (Surface.Apply, BuildMesh, ColliderCook) can be measured uniformly.
    public static class TerrainProfiler
    {
        // -------- Phase identity: the unit at which we attribute time --------

        public enum Phase
        {
            // Main thread
            SurfaceApply,
            ApplyPark,        // Park loop: SetActive(false) presences that left wanted
            ApplyEnsure,      // Ensure loop: revive/create/refresh for chunks in wanted
            ApplyTrim,        // FIFO eviction of warm-cache overflow
            ApplyCount,       // Live-count walk over presences (triangles, etc.)
            BuildMesh,
            MeshAssign,
            ColliderCook,
            DrainResults,
            InstallChunk,
            // Worker thread (timing per "frame" still means "ms of worker work since
            // last main-thread snapshot" — useful for spotting bursts of background work)
            RunPass,
            GenerateChunk,
            GenPrimal,
            GenElevation,
            GenClassify,
            RelaxBorders,
            BuildDual,
            BuildDualGather,    // sub-phase: walk 6 neighbours to build neighborFaces/minNeighborCoord
            BuildDualGenerate,  // sub-phase: the actual GenerateDual polygon math
            DeepCopy,
            Count
        }

        public static readonly bool[] IsMainThreadPhase = BuildMainThreadMask();

        static bool[] BuildMainThreadMask()
        {
            var m = new bool[(int)Phase.Count];
            m[(int)Phase.SurfaceApply]   = true;
            m[(int)Phase.ApplyPark]      = true;
            m[(int)Phase.ApplyEnsure]    = true;
            m[(int)Phase.ApplyTrim]      = true;
            m[(int)Phase.ApplyCount]     = true;
            m[(int)Phase.BuildMesh]      = true;
            m[(int)Phase.MeshAssign]     = true;
            m[(int)Phase.ColliderCook]   = true;
            m[(int)Phase.DrainResults]   = true;
            m[(int)Phase.InstallChunk]   = true;
            return m;
        }

        public static string Name(Phase p) => Names[(int)p];

        static readonly string[] Names =
        {
            "Surface.Apply", "Apply.Park", "Apply.Ensure", "Apply.Trim", "Apply.Count",
            "BuildMesh", "MeshAssign", "ColliderCook",
            "DrainResults", "InstallChunk",
            "RunPass", "GenerateChunk", "Gen.Primal", "Gen.Elevation", "Gen.Classify",
            "RelaxBorders", "BuildDual", "BuildDual.Gather", "BuildDual.Generate", "DeepCopy",
        };

        // -------- Markers (Unity Profiler) --------
        // One ProfilerMarker per Phase, so the Profiler window names match the overlay rows.

        static readonly ProfilerMarker[] markers = BuildMarkers();

        static ProfilerMarker[] BuildMarkers()
        {
            var arr = new ProfilerMarker[(int)Phase.Count];
            for (int i = 0; i < arr.Length; i++)
                arr[i] = new ProfilerMarker("Terrain." + Names[i]);
            return arr;
        }

        // -------- Cumulative ticks per phase --------
        // Monotonically increasing. The overlay snapshots them at frame boundaries and
        // diffs to compute "ms spent in phase X this frame" (or "since last frame snapshot"
        // for worker phases).

        static readonly long[] cumulativeTicks = new long[(int)Phase.Count];

        public static long ReadCumulativeTicks(Phase p) => Interlocked.Read(ref cumulativeTicks[(int)p]);

        public static double TicksToMs(long ticks) => ticks * 1000.0 / Stopwatch.Frequency;

        // The canonical scoping helper. Wraps the corresponding ProfilerMarker AND a
        // Stopwatch-based tick accumulator in one struct, so call sites stay one-liners:
        //     using var _ = TerrainProfiler.Measure(Phase.BuildMesh);
        public struct PerfScope : IDisposable
        {
            int phaseIndex;
            long startTicks;
            bool armed;

            internal PerfScope(int phaseIndex)
            {
                this.phaseIndex = phaseIndex;
                markers[phaseIndex].Begin();
                startTicks = Stopwatch.GetTimestamp();
                armed = true;
            }

            public void Dispose()
            {
                if (!armed) return;
                armed = false;
                long elapsed = Stopwatch.GetTimestamp() - startTicks;
                Interlocked.Add(ref cumulativeTicks[phaseIndex], elapsed);
                markers[phaseIndex].End();
            }
        }

        public static PerfScope Measure(Phase p) => new PerfScope((int)p);

        // -------- Throughput counters (unchanged) --------

        static long _chunksGenerated;
        static long _dualsBuilt;
        static long _relaxPasses;
        static long _deepCopies;
        static long _meshRebuilds;
        static long _colliderAssigns;
        static long _surfaceApplies;
        static long _trianglesBuilt;

        static int _loadedChunks;
        static int _cachedChunks;
        static int _presences;
        static int _presencesCached;
        static int _liveTriangles;

        public static long ChunksGenerated  => Interlocked.Read(ref _chunksGenerated);
        public static long DualsBuilt       => Interlocked.Read(ref _dualsBuilt);
        public static long RelaxPasses      => Interlocked.Read(ref _relaxPasses);
        public static long DeepCopies       => Interlocked.Read(ref _deepCopies);
        public static long MeshRebuilds     => Interlocked.Read(ref _meshRebuilds);
        public static long ColliderAssigns  => Interlocked.Read(ref _colliderAssigns);
        public static long SurfaceApplies   => Interlocked.Read(ref _surfaceApplies);
        public static long TrianglesBuilt   => Interlocked.Read(ref _trianglesBuilt);

        public static int LoadedChunks      { get => Volatile.Read(ref _loadedChunks);   set => Volatile.Write(ref _loadedChunks, value); }
        public static int CachedChunks      { get => Volatile.Read(ref _cachedChunks);   set => Volatile.Write(ref _cachedChunks, value); }
        public static int Presences         { get => Volatile.Read(ref _presences);      set => Volatile.Write(ref _presences, value); }
        public static int PresencesCached   { get => Volatile.Read(ref _presencesCached); set => Volatile.Write(ref _presencesCached, value); }
        public static int LiveTriangles     { get => Volatile.Read(ref _liveTriangles);  set => Volatile.Write(ref _liveTriangles, value); }

        public static void IncChunksGenerated()  => Interlocked.Increment(ref _chunksGenerated);
        public static void IncDualsBuilt()       => Interlocked.Increment(ref _dualsBuilt);
        public static void IncRelaxPasses()      => Interlocked.Increment(ref _relaxPasses);
        public static void IncDeepCopies()       => Interlocked.Increment(ref _deepCopies);
        public static void IncMeshRebuilds()     => Interlocked.Increment(ref _meshRebuilds);
        public static void IncColliderAssigns()  => Interlocked.Increment(ref _colliderAssigns);
        public static void IncSurfaceApplies()   => Interlocked.Increment(ref _surfaceApplies);
        public static void AddTrianglesBuilt(int n) => Interlocked.Add(ref _trianglesBuilt, n);
    }
}
