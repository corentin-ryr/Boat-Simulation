# Terrain Elevation — Alternative Approaches

## Scope

The shipped elevation system (see `TERRAIN_SYSTEM.md` once integrated) is intentionally the
simplest pointwise field that produces island-dotted oceans:

```
height(x, z) = smoothstep(t, t + ε, fBm(x, z)) · profile(fBm(x, z))
```

— a single OpenSimplex2 fBm field, soft-thresholded so most of the world is exactly zero
(flat ocean), with islands rising above. This document catalogs the techniques that were
*considered and deferred*. Each is pointwise (no simulation, no state), so any of them can
be dropped in later without touching the streaming/meshing/collider layers — only the
implementation of `Elevation.Sample(x, z)` changes.

The flat-chunk optimization (`PrimalChunk.IsFlat`) is preserved by every option below, as
long as the technique reduces to exactly zero below threshold.

---

## 1. Island Mask × Elevation Profile (two-layer)

**Idea.** Split "where can land exist?" from "what does the land look like?":

```
mask     = fBm_low (x, z)              // few octaves, large features
profile  = fBm_high(x, z)              // many octaves, fine detail
landness = smoothstep(t, t+ε, mask)    // 0 in ocean, 1 well inland
height   = landness · profile_curve(profile)
```

**Why upgrade.** With a single field, islands form a uniform Perlin-spatter pattern across
the entire map. Splitting the mask out lets you tune **island clustering** (large quiet
zones vs. dense archipelagos) independently from **island shape** (smooth domes vs.
craggy plateaus). Two independent knobs.

**Cost.** ~2× noise calls per sample. Still cheap.

**Variants.**
- Use **distance-from-mask-edge** instead of mask amplitude as the landness term — gives
  more uniform island heights regardless of how "deep" inland you are.
- Multiply by a coarse **continental mask** at world scale so the whole map has macro
  regions of "ocean-heavy" vs. "land-heavy."

**References.**
- Sebastian Lague — *Procedural Planet Generation* series:
  https://www.youtube.com/playlist?list=PLFt_AvWsXl0cONs3T0By4puYy6GM22ko8
- Catlike Coding — *Noise* tutorial (still the best step-by-step for layered noise in Unity):
  https://catlikecoding.com/unity/tutorials/noise/

---

## 2. Mask + Profile + Ridged Mountains (three-layer)

**Idea.** Add a third layer on top of (1) that introduces sharp ridges only where the
landmass is thick:

```
ridged   = 1 - |fBm_ridged(x, z)|         // ridges where noise crosses zero
mountains = smoothstep(m_lo, m_hi, mask) · ridged
height   = (landness · profile_curve(profile)) + α · mountains
```

**Why upgrade.** Real archipelagos have a wide spectrum: tiny flat sandbars, rolling
medium islands, and a few mountainous ones. Ridged multifractal noise gives you the
crisp peaks that fBm can't — fBm's spectrum is too smooth at the top, mountains end up
rounded.

**Cost.** ~3–4× noise calls per sample, plus a couple of smoothsteps. Still cheap; this
is well within budget for per-vertex evaluation at dual-build time.

**Tuning trap.** Ridged noise on its own looks alien — like spiderwebs of ridges over
*everything*. Always gate it by the mask amplitude so it only appears where land is
already established.

**References.**
- Inigo Quilez — *Fractals: Ridged multifractal*: https://iquilezles.org/articles/morenoise/
- Original paper: Musgrave, *"Methods for Realistic Landscape Imaging"* (1993).
  Chapter 16 of *Texturing & Modeling: A Procedural Approach* (Ebert et al., 3rd ed.)
  is the canonical writeup.

---

## 3. Domain Warping

**Idea.** Don't sample noise at `(x, z)`; sample it at `(x + q(x,z), z + r(x,z))`, where
`q` and `r` are themselves noise fields:

```
qx = fBm(x + 0,   z + 0)
qz = fBm(x + 5.2, z + 1.3)
height = fBm(x + α·qx, z + α·qz)
```

**Why upgrade.** Raw fBm produces "fingerprint" blobs with self-similar swirls at every
scale — pleasant but obviously synthetic. One pass of domain warping produces curvy,
flowing coastlines reminiscent of Norwegian fjords or Greek islands. It's the single
biggest **aesthetic win per CPU cycle** in procedural terrain.

**Cost.** 2–3 extra noise evaluations per sample.

**Stacking.** Iterated warping (`noise(noise(noise(...)))`) is the iconic Inigo Quilez
trick. Each level adds another layer of organic flow at the cost of N more noise calls.
Two levels is usually the sweet spot.

**References.**
- Inigo Quilez — *Domain warping* (the canonical writeup, with stunning examples):
  https://iquilezles.org/articles/warp/
- Shadertoy demonstrations are linked from that article.

---

## 4. Worley / Voronoi-Based Features

**Idea.** Use cellular noise (distance-to-nearest-feature-point) as a secondary layer:

```
F1     = distance to nearest random point
F2     = distance to 2nd nearest
ridges = F2 - F1     // small where two cells nearly meet → cracks
cells  = F1          // smooth domes per cell
```

**Why upgrade.** Adds discrete, structured features that pure-noise approaches miss:
- **River-like cracks** along Voronoi edges (use `F2 - F1 < ε`).
- **Crystalline rock outcrops** (use `F1` modulated by mask).
- **Discrete biome regions** that map naturally to your hex tiles.

**Cost.** A bit more than fBm — distance to N nearest points requires querying a 3×3 cell
neighborhood per sample, but each cell only needs one or two points.

**Synergy with your hex grid.** Your primal grid is already cell-based. You can use the
primal cell index itself as the Voronoi seed and skip the noise lookup entirely — every
hex naturally becomes a "biome cell," which dovetails with the planned "1 hex = 1 tile
of meaning" gameplay design.

**References.**
- Worley's original paper: *"A Cellular Texture Basis Function"* (1996).
- Inigo Quilez — *Voronoi distances*: https://iquilezles.org/articles/voronoilines/
- KdotJPG's OpenSimplex2 repo ships a cellular variant:
  https://github.com/KdotJPG/OpenSimplex2

---

## 5. Continental Falloff

**Idea.** Multiply the mask by a coarse falloff function based on world position to
constrain land to a specific macro region:

```
continental = falloff(distance_to_world_center)
mask = continental · fBm_low(x, z)
```

`falloff` could be a radial smoothstep (single round continent), a hand-painted texture
sampled from a low-res map (designer-controlled coastline), or another low-frequency
noise at planetary scale (multi-continent worlds).

**Why use.** Lets you guarantee "the player spawn area is always ocean with islands
within 5 km" without leaving the world's geometry entirely up to noise. Especially
useful when you want a bounded play area instead of an infinite procedural map.

**Cost.** Effectively free (one extra sample, no octaves).

**References.**
- Sebastian Lague — *Procedural Coastlines* (uses falloff maps explicitly):
  https://www.youtube.com/watch?v=COZHVUPp_qg

---

## 6. Out of Scope — Why These Are Not Pointwise

Listed for completeness; rejected because they require global state or simulation:

- **Hydraulic / thermal erosion.** Drops virtual rain on a heightmap and simulates flow
  and deposition. Produces stunningly realistic terrain (the GTA / Horizon look) but
  needs the entire heightmap at once — incompatible with streaming a chunk on demand.
  Could in principle be baked offline per-chunk with halo overlap, but the seams are
  hard.
- **Tectonic plate simulation.** Discretizes the world into plates and simulates their
  collision over time. Produces realistic mountain ranges along plate boundaries.
  Inherently global.
- **Diamond-square / midpoint displacement.** Recursive subdivision of a coarse grid.
  Was the classic 90s-era approach; predates noise. Not pointwise — each sample depends
  on the recursion path.
- **Wave Function Collapse for tiles.** Produces structured tile-based worlds. Not a
  continuous height field; would replace, not augment, the noise approach.
- **Neural / GAN-based heightmaps.** Train a network on real-world DEMs, sample heights
  via inference. Pointwise in theory, but the inference cost per sample (and the
  inability to back-port the result into the chunk pipeline) makes it impractical for
  per-frame evaluation. Could be used to **bake** a heightmap that's then sampled
  pointwise.

---

## General References

- *Texturing & Modeling: A Procedural Approach* — Ebert, Musgrave, Peachey, Perlin,
  Worley (3rd ed., 2002). The foundational text. Chapters on noise basis functions
  (Perlin/Worley) and procedural terrain (Musgrave's multifractals) are still the
  reference.
- Inigo Quilez — articles index: https://iquilezles.org/articles/
- KdotJPG — OpenSimplex2 reference implementations (C#, public domain):
  https://github.com/KdotJPG/OpenSimplex2
- Red Blob Games — *Making maps with noise functions*:
  https://www.redblobgames.com/maps/terrain-from-noise/
- Sebastian Lague — *Coding Adventures* (terrain, planets, erosion) on YouTube.
- Catlike Coding — *Noise* tutorial series:
  https://catlikecoding.com/unity/tutorials/noise/
