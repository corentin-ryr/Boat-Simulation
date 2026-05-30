using System.Text;
using UnityEngine;

namespace TerrainGrid
{
    // On-screen overlay for diagnosing terrain performance, drawn in both the Game view
    // (OnGUI) and any open Scene view (SceneView.duringSceneGui, editor only).
    //
    // Three sections:
    //
    //   1. Live counts — chunks loaded/cached/presented, total live triangles, frame time
    //      (smoothed). Tells you "is the working set the size I expect?"
    //
    //   2. Peak frame ms by phase — the highest single-frame cost of every instrumented
    //      phase observed in the last `sampleInterval` seconds. THIS is what you watch
    //      when investigating camera-pan slowdowns: a one-frame spike pops out here even
    //      if the average is fine. Main-thread phases (Surface.Apply, BuildMesh, etc.) are
    //      directly responsible for hitches; worker phases (BuildDual, GenerateChunk) show
    //      "ms of background work since the last frame", which is useful for spotting
    //      bursts of generation pressure but doesn't directly stall the main thread.
    //
    //   3. Per-second rates — throughput counters. Tells you "how hard is the pipeline
    //      churning over a moving average?"
    //
    // When a frame's main-thread terrain time exceeds `spikeThresholdMs`, the overlay
    // also logs a one-line breakdown to the console so you can scroll back and see which
    // phases dominated the worst frames without staring at the overlay constantly.
    //
    // Use together with the Unity Profiler: the overlay tells you WHICH phase is the
    // problem; the Profiler tells you WHICH CALL inside that phase was the slow one.
    public class TerrainStatsOverlay : MonoBehaviour
    {
        [Header("Display")]
        [Tooltip("Master toggle. Turn off in builds where you don't want the overlay.")]
        public bool show = true;

        [Tooltip("Also draw inside the editor's Scene view window (editor only).")]
        public bool showInSceneView = true;

        [Tooltip("Seconds between rate samples and peak-window resets. Peak phase ms is the worst single-frame cost observed in this window.")]
        [Min(0.1f)] public float sampleInterval = 0.5f;

        [Tooltip("Overlay corner anchor.")]
        public Corner anchor = Corner.TopRight;

        [Tooltip("Pixel offset from the screen edge.")]
        public Vector2 padding = new Vector2(12f, 12f);

        [Tooltip("Background alpha behind the text (0 = no background).")]
        [Range(0f, 1f)] public float backgroundAlpha = 0.5f;

        [Tooltip("Font size for the overlay text. Bump up on retina/high-DPI displays.")]
        [Range(8, 32)] public int fontSize = 13;

        [Header("Spike Logging")]
        [Tooltip("Log a per-phase breakdown to the console whenever a frame's main-thread terrain time exceeds this many milliseconds.")]
        [Min(0f)] public float spikeThresholdMs = 8f;

        [Tooltip("Maximum spike logs per second so a sustained slowdown doesn't flood the console.")]
        [Min(0f)] public float maxSpikeLogsPerSecond = 4f;

        public enum Corner { TopLeft, TopRight, BottomLeft, BottomRight }

        // -------- internal state --------

        const int PhaseCount = (int)TerrainProfiler.Phase.Count;

        // Last cumulative tick read per phase (for frame deltas).
        readonly long[] prevTicks = new long[PhaseCount];

        // This-frame ms per phase, recomputed every LateUpdate.
        readonly double[] thisFrameMs = new double[PhaseCount];

        // Peak frame ms per phase, updated every frame, reset every sampleInterval.
        readonly double[] peakFrameMs = new double[PhaseCount];
        // Display copy of the peaks, copied at the end of each sampleInterval window.
        readonly double[] displayedPeakMs = new double[PhaseCount];

        // Smoothed frame time, plus peak frame time in the current window.
        float frameMsSmoothed;
        float peakFrameMsTotal;
        float displayedPeakFrameMsTotal;

        // Throughput rates.
        float r_generated, r_dualsBuilt, r_relax, r_deepCopies;
        float r_meshRebuilds, r_colliderAssigns, r_surfaceApplies, r_trianglesBuilt;
        long s_generated, s_dualsBuilt, s_relax, s_deepCopies;
        long s_meshRebuilds, s_colliderAssigns, s_surfaceApplies, s_trianglesBuilt;

        float lastSampleAt;
        float spikeLogBudget;
        float lastSpikeLogTime;

        GUIStyle textStyle;
        GUIStyle boxStyle;
        Texture2D bgTex;

        const float OverlayWidth  = 320f;
        const float OverlayHeight = 460f;

        // Pre-built phase order for display (main thread phases first, then worker).
        static readonly TerrainProfiler.Phase[] DisplayOrder =
        {
            TerrainProfiler.Phase.SurfaceApply,
            TerrainProfiler.Phase.ApplyPark,
            TerrainProfiler.Phase.ApplyEnsure,
            TerrainProfiler.Phase.ApplyTrim,
            TerrainProfiler.Phase.ApplyCount,
            TerrainProfiler.Phase.BuildMesh,
            TerrainProfiler.Phase.MeshAssign,
            TerrainProfiler.Phase.ColliderCook,
            TerrainProfiler.Phase.DrainResults,
            TerrainProfiler.Phase.InstallChunk,
            TerrainProfiler.Phase.RunPass,
            TerrainProfiler.Phase.GenerateChunk,
            TerrainProfiler.Phase.GenPrimal,
            TerrainProfiler.Phase.GenElevation,
            TerrainProfiler.Phase.GenClassify,
            TerrainProfiler.Phase.RelaxBorders,
            TerrainProfiler.Phase.BuildDual,
            TerrainProfiler.Phase.BuildDualGather,
            TerrainProfiler.Phase.BuildDualGenerate,
            TerrainProfiler.Phase.DeepCopy,
        };

        void OnEnable()
        {
            // Prime previous-tick snapshot so the first frame's delta is sensible.
            for (int i = 0; i < PhaseCount; i++)
                prevTicks[i] = TerrainProfiler.ReadCumulativeTicks((TerrainProfiler.Phase)i);

#if UNITY_EDITOR
            UnityEditor.SceneView.duringSceneGui += OnSceneViewGUI;
#endif
        }

        void OnDisable()
        {
#if UNITY_EDITOR
            UnityEditor.SceneView.duringSceneGui -= OnSceneViewGUI;
#endif
        }

        // LateUpdate runs after every other component has finished its frame, including
        // ChunkSurface.Apply (driven from ChunkManager.LateUpdate). So the cumulative tick
        // delta we read here is exactly "ms spent in each phase during this frame".
        void LateUpdate()
        {
            // 1) Frame time (smoothed for the live readout, raw for spike detection).
            float dtMs = Time.unscaledDeltaTime * 1000f;
            frameMsSmoothed = Mathf.Lerp(frameMsSmoothed, dtMs, 0.1f);
            if (dtMs > peakFrameMsTotal) peakFrameMsTotal = dtMs;

            // 2) Per-phase delta this frame. Updates thisFrameMs (cleared first) and
            //    folds each phase's value into the rolling peak.
            for (int i = 0; i < PhaseCount; i++)
            {
                long now = TerrainProfiler.ReadCumulativeTicks((TerrainProfiler.Phase)i);
                long delta = now - prevTicks[i];
                prevTicks[i] = now;
                double ms = TerrainProfiler.TicksToMs(delta);
                thisFrameMs[i] = ms;
                if (ms > peakFrameMs[i]) peakFrameMs[i] = ms;
            }

            // Total main-thread terrain ms this frame = sum of the umbrella phases only
            // (SurfaceApply already contains BuildMesh/MeshAssign/ColliderCook; DrainResults
            // already contains InstallChunk). Adding leaves would double-count.
            double mainTotalThisFrame =
                  thisFrameMs[(int)TerrainProfiler.Phase.SurfaceApply]
                + thisFrameMs[(int)TerrainProfiler.Phase.DrainResults];

            // 3) Spike detection: if the main-thread terrain spend this frame is over the
            //    threshold, log a one-line breakdown (rate-limited).
            if (spikeThresholdMs > 0f && mainTotalThisFrame >= spikeThresholdMs)
                TryLogSpike(dtMs, mainTotalThisFrame);

            // 4) Throughput sampling on the slow timer.
            float now2 = Time.unscaledTime;
            float window = now2 - lastSampleAt;
            if (window >= sampleInterval)
            {
                lastSampleAt = now2;

                long g  = TerrainProfiler.ChunksGenerated;
                long db = TerrainProfiler.DualsBuilt;
                long rp = TerrainProfiler.RelaxPasses;
                long dc = TerrainProfiler.DeepCopies;
                long mr = TerrainProfiler.MeshRebuilds;
                long ca = TerrainProfiler.ColliderAssigns;
                long sa = TerrainProfiler.SurfaceApplies;
                long tb = TerrainProfiler.TrianglesBuilt;
                float invDt = 1f / window;
                r_generated       = (g  - s_generated)       * invDt;
                r_dualsBuilt      = (db - s_dualsBuilt)      * invDt;
                r_relax           = (rp - s_relax)           * invDt;
                r_deepCopies      = (dc - s_deepCopies)      * invDt;
                r_meshRebuilds    = (mr - s_meshRebuilds)    * invDt;
                r_colliderAssigns = (ca - s_colliderAssigns) * invDt;
                r_surfaceApplies  = (sa - s_surfaceApplies)  * invDt;
                r_trianglesBuilt  = (tb - s_trianglesBuilt)  * invDt;
                s_generated = g; s_dualsBuilt = db; s_relax = rp; s_deepCopies = dc;
                s_meshRebuilds = mr; s_colliderAssigns = ca; s_surfaceApplies = sa;
                s_trianglesBuilt = tb;

                // Snapshot the window's peaks to the display arrays, then reset.
                for (int i = 0; i < PhaseCount; i++)
                {
                    displayedPeakMs[i] = peakFrameMs[i];
                    peakFrameMs[i] = 0;
                }
                displayedPeakFrameMsTotal = peakFrameMsTotal;
                peakFrameMsTotal = 0;
            }

            // 5) Spike-log rate limiter recharge.
            if (maxSpikeLogsPerSecond > 0f)
                spikeLogBudget = Mathf.Min(spikeLogBudget + Time.unscaledDeltaTime * maxSpikeLogsPerSecond,
                                            maxSpikeLogsPerSecond);
        }

        void TryLogSpike(float frameMs, double mainTotalMs)
        {
            if (maxSpikeLogsPerSecond > 0f && spikeLogBudget < 1f) return;
            spikeLogBudget -= 1f;

            var sb = new StringBuilder(256);
            sb.Append("[Terrain spike] frame ").Append(frameMs.ToString("0.0"))
              .Append("ms  main-thread terrain ").Append(mainTotalMs.ToString("0.0")).Append("ms — ");

            bool first = true;
            for (int i = 0; i < PhaseCount; i++)
            {
                if (!TerrainProfiler.IsMainThreadPhase[i]) continue;
                double ms = thisFrameMs[i];
                if (ms < 0.1) continue;
                if (!first) sb.Append(", ");
                first = false;
                sb.Append(TerrainProfiler.Name((TerrainProfiler.Phase)i))
                  .Append(' ').Append(ms.ToString("0.0")).Append("ms");
            }

            Debug.Log(sb.ToString());
            lastSpikeLogTime = Time.unscaledTime;
        }

        // -------- rendering --------

        void OnGUI()
        {
            if (!show) return;
            DrawOverlay(Screen.width, Screen.height);
        }

#if UNITY_EDITOR
        void OnSceneViewGUI(UnityEditor.SceneView sv)
        {
            if (!show || !showInSceneView) return;
            UnityEditor.Handles.BeginGUI();
            DrawOverlay((int)sv.position.width, (int)sv.position.height);
            UnityEditor.Handles.EndGUI();
        }
#endif

        void DrawOverlay(int screenW, int screenH)
        {
            EnsureStyles();

            float x = anchor == Corner.TopLeft || anchor == Corner.BottomLeft
                ? padding.x
                : screenW - OverlayWidth - padding.x;
            float y = anchor == Corner.TopLeft || anchor == Corner.TopRight
                ? padding.y
                : screenH - OverlayHeight - padding.y;

            var rect = new Rect(x, y, OverlayWidth, OverlayHeight);

            if (backgroundAlpha > 0f) GUI.Box(rect, GUIContent.none, boxStyle);

            GUILayout.BeginArea(new Rect(x + 8, y + 6, OverlayWidth - 16, OverlayHeight - 12));

            GUILayout.Label("<b>Terrain</b>", textStyle);
            GUILayout.Label($"frame: {frameMsSmoothed,5:0.0} ms ({1000f / Mathf.Max(frameMsSmoothed, 0.01f),4:0} fps)  peak: {displayedPeakFrameMsTotal:0.0}ms", textStyle);
            GUILayout.Space(2f);
            GUILayout.Label($"loaded: {TerrainProfiler.LoadedChunks,4}   cached: {TerrainProfiler.CachedChunks,3}", textStyle);
            GUILayout.Label($"presented: {TerrainProfiler.Presences,4}   warm: {TerrainProfiler.PresencesCached,3}", textStyle);
            GUILayout.Label($"tris: {FormatThousands(TerrainProfiler.LiveTriangles)}", textStyle);

            GUILayout.Space(4f);
            GUILayout.Label("<b>peak frame ms (last interval)</b>", textStyle);
            GUILayout.Label("<i>main thread:</i>", textStyle);
            DrawPhasePeak(TerrainProfiler.Phase.SurfaceApply);
            DrawPhasePeak(TerrainProfiler.Phase.BuildMesh,    indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.MeshAssign,   indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.ColliderCook, indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.DrainResults);
            DrawPhasePeak(TerrainProfiler.Phase.InstallChunk, indent: true);
            GUILayout.Label("<i>worker:</i>", textStyle);
            DrawPhasePeak(TerrainProfiler.Phase.RunPass);
            DrawPhasePeak(TerrainProfiler.Phase.GenerateChunk, indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.BuildDual,    indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.BuildDualGather,   indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.BuildDualGenerate, indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.RelaxBorders, indent: true);
            DrawPhasePeak(TerrainProfiler.Phase.DeepCopy,     indent: true);

            GUILayout.Space(4f);
            GUILayout.Label("<b>per second</b>", textStyle);
            GUILayout.Label($"  gen {r_generated,5:0.0}   dual {r_dualsBuilt,5:0.0}   relax {r_relax,5:0.0}", textStyle);
            GUILayout.Label($"  mesh {r_meshRebuilds,5:0.0}  collider {r_colliderAssigns,5:0.0}  copy {r_deepCopies,4:0.0}", textStyle);

            GUILayout.EndArea();
        }

        void DrawPhasePeak(TerrainProfiler.Phase p, bool indent = false)
        {
            double ms = displayedPeakMs[(int)p];
            string label = (indent ? "  " : "") + TerrainProfiler.Name(p);
            string color = ms >= spikeThresholdMs ? "#ff8080"
                         : ms >= spikeThresholdMs * 0.5f ? "#ffe080"
                         : "#a0e0a0";
            GUILayout.Label($"<color={color}>{label,-20} {ms,6:0.0} ms</color>", textStyle);
        }

        void EnsureStyles()
        {
            if (textStyle == null || textStyle.fontSize != fontSize)
            {
                textStyle = new GUIStyle(GUI.skin.label)
                {
                    fontSize = fontSize,
                    richText = true,
                    normal = { textColor = Color.white },
                };
            }
            if (boxStyle == null || bgTex == null)
            {
                bgTex = new Texture2D(1, 1) { hideFlags = HideFlags.HideAndDontSave };
                bgTex.SetPixel(0, 0, new Color(0f, 0f, 0f, backgroundAlpha));
                bgTex.Apply();
                boxStyle = new GUIStyle(GUI.skin.box) { normal = { background = bgTex } };
            }
        }

        void OnDestroy()
        {
            if (bgTex != null) Destroy(bgTex);
        }

        static string FormatThousands(int n)
        {
            if (n >= 1_000_000) return $"{n / 1_000_000f:0.0}M";
            if (n >= 1_000)     return $"{n / 1_000f:0.0}K";
            return n.ToString();
        }
    }
}
