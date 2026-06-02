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
    // The field is bit-exactly `SeabedY` (= -seabedDepth) below the threshold, so chunks that
    // are entirely under it can be flagged IsFlat and short-circuit dual generation, mesh build,
    // and per-vertex storage. The seabed sits below the (Y=0) ocean surface so the FFT water
    // surface visibly covers it from above. See ELEVATION_ALTERNATIVES.md for alternatives.
    public static class Elevation
    {
        // Assigned once at startup (typically from ChunkManager.Start). Worker threads then
        // read it without locking — config is intended to be effectively immutable post-init.
        public static ElevationConfig Config;

        // World Y of the flat seabed (= -seabedDepth). Single source of truth used by both the
        // sampler and ChunkSurface's flat-tile placement so they cannot drift apart.
        public static float SeabedY => -(Config?.seabedDepth ?? 0f);

        // Height at world position (x, z).
        //
        // Base shape: linearly map noise to height as `(noise - threshold) × heightScale`, then
        // clamp to the seabed level — `h = max((noise - threshold) × heightScale, -seabedDepth)`.
        // The shoreline (h crosses 0) is in the interior of the linear region, so it inherits
        // fBm's C∞ smoothness automatically — no shoreline smoothing needed.
        //
        // The remaining C¹ kink is underwater, where the linear slope meets the flat seabed
        // clamp. Optional: set `beachBand > 0` to smooth that crease via a smoothstep over
        // `beachBand` noise units worth of slope. C² on the seabed side, C¹ on the slope side.
        // Set `beachBand = 0` to keep the sharp `max` and rely on water hiding the kink.
        //
        // Samples that fall in the pure-seabed region return `-seabedDepth` bit-exactly, so
        // the IsFlat fast path stays bit-exact and chunks fully below the slope keep using the
        // shared flat tile.
        public static float Sample(float x, float z)
        {
            ElevationConfig cfg = Config;
            if (cfg == null) return 0f;
            float seabed = -cfg.seabedDepth;
            if (cfg.octaves <= 0 || cfg.heightScale <= 0f) return seabed;

            float n = FBm(x, z, cfg);
            float h = (n - cfg.threshold) * cfg.heightScale;

            // Sharp clamp: pure max(h, seabed). Cheapest path, C¹ kink at the seabed boundary.
            if (cfg.beachBand <= 0f) return h <= seabed ? seabed : h;

            // Smoothed clamp: lerp from seabed up to the linear slope over `bandHeight` worth
            // of vertical distance. The smoothstep has s(0)=s'(0)=0 and s(1)=1, s'(1)=0, so:
            //   - At the seabed boundary, value and slope both come in at zero → C¹ join with
            //     the flat clamp below. The product `s × (h-seabed)` also kills the second
            //     derivative there → C² on the seabed side.
            //   - At the upper edge of the band, s saturates at 1 → seamless C¹ join with the
            //     linear slope above.
            float bandHeight = cfg.beachBand * cfg.heightScale;
            float t = (h - seabed) / bandHeight;
            if (t <= 0f) return seabed;
            if (t >= 1f) return h;
            return Mathf.Lerp(seabed, h, SmoothStep01(t));
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

        [Tooltip("Optional underwater smoothing of the seabed→slope crease, in noise units. " +
                 "0 = sharp clamp (cheapest; kink is invisible under water anyway). >0 = " +
                 "smoothstep blend over this much slope, giving a C¹ rounded contour where " +
                 "the linear bathymetry meets the flat seabed. The shoreline self-smooths " +
                 "in either case (fBm is C∞ at the waterline by construction).")]
        [Min(0f)] public float beachBand = 0f;

        [Tooltip("World units of elevation per unit of (noise - threshold). Bigger = taller " +
                 "islands. With threshold=0.15 and noise in [-1,1], an island peak reaches " +
                 "about 0.85 × heightScale.")]
        [Min(0f)] public float heightScale = 20f;

        [Tooltip("How far the flat ocean bed sits below water level (Y=0), in world units. " +
                 "Increase to make the seabed visibly deeper through the water; set 0 to " +
                 "match the old coplanar behaviour. The beach band smoothly lerps from this " +
                 "depth up to the land profile so coastlines don't tear at the threshold.")]
        [Min(0f)] public float seabedDepth = 5f;

        [Tooltip("Independent noise seed. Change to get a different world layout without " +
                 "changing the chunk RNG seed (which controls triangle merge patterns).")]
        public long seed = 1337L;
    }
}
