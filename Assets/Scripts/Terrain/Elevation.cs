using System;
using TerrainGrid.Noise;
using UnityEngine;

namespace TerrainGrid
{
    // Pointwise terrain height field. Pure function of (x, z) — no state, thread-safe — so the
    // same call can be made from the worker thread (dual mesh build, IsFlat classification),
    // the main thread (collider rebuild), and gameplay queries (NPC ground snap, building
    // placement) and every caller gets the same answer.
    //
    // The field is exactly zero below the threshold, so chunks that are entirely under it can
    // be flagged IsFlat and short-circuit dual generation, mesh build, and per-vertex storage.
    // See ELEVATION_ALTERNATIVES.md for the alternatives that were considered.
    public static class Elevation
    {
        // Assigned once at startup (typically from ChunkManager.Start). Worker threads then
        // read it without locking — config is intended to be effectively immutable post-init.
        public static ElevationConfig Config;

        // Height at world position (x, z). Returns exactly 0 for any sample whose noise value
        // is below the configured threshold — that's what makes the flat-ocean fast path
        // bit-exact rather than "approximately zero."
        public static float Sample(float x, float z)
        {
            ElevationConfig cfg = Config;
            if (cfg == null || cfg.octaves <= 0 || cfg.heightScale <= 0f) return 0f;

            float n = FBm(x, z, cfg);
            float above = n - cfg.threshold;
            if (above <= 0f) return 0f;

            // Smooth ramp over `beachBand` for soft shorelines. Above the band the ramp
            // saturates and elevation is purely the profile term.
            float gate = cfg.beachBand > 0f && above < cfg.beachBand
                ? SmoothStep01(above / cfg.beachBand)
                : 1f;

            return gate * above * cfg.heightScale;
        }

        // Fractal Brownian Motion — sum of `octaves` noise samples at successively higher
        // frequencies and lower amplitudes. Normalized so the output stays in roughly [-1, 1]
        // regardless of octave count, making the threshold value meaningful across configs.
        static float FBm(float x, float z, ElevationConfig cfg)
        {
            float amp = 1f, freq = cfg.baseFrequency;
            float sum = 0f, norm = 0f;
            long seed = cfg.seed;
            for (int o = 0; o < cfg.octaves; o++)
            {
                sum += amp * OpenSimplex2.Noise2(seed + o, x * freq, z * freq);
                norm += amp;
                amp *= cfg.gain;
                freq *= cfg.lacunarity;
            }
            return norm > 0f ? sum / norm : 0f;
        }

        // Standard cubic smoothstep on [0, 1]. Returns 0 at t<=0, 1 at t>=1, smooth in between.
        static float SmoothStep01(float t)
        {
            if (t <= 0f) return 0f;
            if (t >= 1f) return 1f;
            return t * t * (3f - 2f * t);
        }
    }

    // Tuning knobs for the elevation field. Inlined on ChunkManager (not a ScriptableObject)
    // for v1 simplicity — promote to an asset later if multiple scenes need different presets.
    [Serializable]
    public class ElevationConfig
    {
        [Tooltip("Number of fBm octaves. ~4–6 is the conventional range. Each octave costs one " +
                 "noise call per sample (per primal vertex, computed once at chunk generation).")]
        [Min(0)] public int octaves = 5;

        [Tooltip("Base frequency of the lowest octave, in cycles per world unit. Smaller = " +
                 "larger islands. 0.01 means islands ~100 units across at the macro scale.")]
        public float baseFrequency = 0.01f;

        [Tooltip("Frequency multiplier per octave. Conventional 2.0 doubles the frequency each " +
                 "step (and halves the wavelength).")]
        public float lacunarity = 2.0f;

        [Tooltip("Amplitude multiplier per octave. Conventional 0.5 halves the amplitude each " +
                 "step, so higher octaves contribute fine detail without dominating the shape.")]
        public float gain = 0.5f;

        [Tooltip("Noise value below which terrain is flat ocean (exactly Y=0). Raising this " +
                 "shrinks the land fraction; lowering grows it. Range [-1, 1].")]
        [Range(-1f, 1f)] public float threshold = 0.15f;

        [Tooltip("Soft transition width on the land side of the threshold (in noise units). " +
                 "Set 0 for hard cliffs at every shoreline; small positive values give gentle " +
                 "beaches. Tune in concert with heightScale to control the beach slope.")]
        [Min(0f)] public float beachBand = 0.05f;

        [Tooltip("World units of elevation per unit of (noise - threshold). Bigger = taller " +
                 "islands. With threshold=0.15 and noise in [-1,1], an island peak reaches " +
                 "about 0.85 × heightScale.")]
        [Min(0f)] public float heightScale = 20f;

        [Tooltip("Independent noise seed. Change to get a different world layout without " +
                 "changing the chunk RNG seed (which controls triangle merge patterns).")]
        public long seed = 1337L;
    }
}
