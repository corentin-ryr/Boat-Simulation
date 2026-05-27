# Boat Simulation — Project Architecture

Unity 6 (6000.3) project using URP. A boat physics simulation combining spectral ocean generation, GPU-accelerated buoyancy, and clipmap LOD rendering.

## System Overview

```
WavesGenerator (FFT, GPU)
    └─> WaterManager.GetWaterHeight()
            └─> Floater / PressureFloater
                    ├─> FloaterHelper (waterline trace)
                    ├─> MeshHelper (geometry)
                    └─> Rigidbody.AddForce()

OceanGeometry (LOD mesh) ──> Ocean Shader ──> PlanarReflections
                                           └─> DepthNormalsFeature
```

## Water Surface

`IWater` interface — four interchangeable implementations selected per scene:
- `FlatWater` — static Y=0 plane
- `SinusoidalWater` — simple `sin(x + time)` waves
- `WaterGerstner` — superposition of 4 Gerstner waves
- `WavesGenerator` — full FFT ocean (JONSWAP/Phillips spectrum, 3 GPU cascade levels, async readback for physics)

`WaterManager` is the single access point for all height queries.

## Buoyancy Physics

Two approaches (swap per scene/object):

**`Floater.cs`** — volume-based (physically accurate)
- 5×5×5 spatial grid over boat mesh
- GPU compute shader finds boat/water triangle overlaps
- `FloaterHelper.ComputeIntersectingLine()` traces the exact waterline (flood-fill)
- Integrates pressure over submerged volume → `Rigidbody.AddForce()`
- `RealTimeVector3DSmoother` reduces jitter

**`PressureFloater.cs`** — vertex pressure-based (simpler, more stable)
- Depth per vertex below waterline
- 3 geometric split-triangle cases at waterline
- Viscous linear drag + angular drag (rho=1000, mu=1e-3)

## Ocean Rendering

- **`OceanGeometry.cs`** — clipmap LOD: central mesh + up to 8 concentric rings, grid-snapped positions, shader keywords `CLOSE`/`MID` by camera height
- **`PlanarReflections.cs`** — mirror-camera reflections with oblique clip plane → `_PlanarReflectionTexture`
- **`DepthNormalsFeature.cs`** — custom URP renderer feature; re-renders geometry with normals override material → `_CameraDepthNormalsTexture` for outline shader

## Helper Layer

| File | Role |
|---|---|
| `MeshHelper.cs` | Neighbor graphs, vertex welding, spatial grids, volume computation |
| `FloaterHelper.cs` | Waterline intersection tracing (edge/plane intersection cases) |
| `DataStructures.cs` | `Triangle`, `Vertex`, `Cell`, `MovingAverage` |
| `ComputeHelper.cs` | Compute shader abstraction: thread groups, buffers, reflection-based param binding |
| `FastFourierTransform.cs` | 2D GPU FFT with ping-pong buffers and twiddle factor precomputation |
| `Isocahedron.cs` | GPU-subdivided procedural sphere for test floating objects |
| `Geometry.cs` | `Edge` struct for icosahedron mesh topology |
| `DebugHelper.cs` | Gizmo drawing for forces, bounds, wireframe meshes |

## Terrain Generation (Townscaper-style primal/dual)

`Assets/Scripts/Terrain/` — an irregular quad-grid terrain built with a deliberate
**primal/dual split**. See `TERRAIN_SYSTEM.md` for the full pipeline.

**Primal grid** (`PolygonGridGenerator.GeneratePrimal`) — the *logical/data* grid.
Townscaper construction: triangle lattice (`GenHexShape`) → random triangle merge →
ortho subdivide into quads (`SubdividePolygon`) → relax (borders pinned). Cheap to
generate, so it can be produced broadly — including non-visible parts of the map —
for occupancy, adjacency, gameplay, and streaming logic.

**Dual grid** (`GenerateDualGrid`) — the *rendered surface + transitions* grid.
Standard mesh dual (primal vertex → dual cell, primal face-center → dual vertex), so
dual cells sit on the junctions where primal cells meet (natural home for
corner/transition meshing). More expensive, so compute it **only for the near-camera
region**. The dual is what gets rendered.

Design intent: primal answers "what exists?" (cheap, wide), dual answers "how does it
render/connect?" (costly, only where visible).

**Chunking** (`ChunkCoord`) — each chunk is a flat-top regular hexagon region of the
triangle lattice, circumradius `chunkGridSize · hexRadius`; chunks tile as a honeycomb
(`WorldCenter` uses `size = chunkGridSize · hexRadius`). No ghost ring — border vertices
coincide exactly across neighbors. Cross-chunk matching (seam welding in `RelaxBorders`,
dual cell completion + ownership in `GenerateDual`) keys on a deterministic integer
`LatticeKey(position)`, not rounded floats, so coincident copies always match.

**Streaming / threading** — decoupled layers, each in its own file:
- `TerrainModel` — the authoritative, Unity-free data layer: owns loaded `PrimalChunk`s,
  generation, joint border relaxation, and the chunk store.
- `IChunkStore` / `MemoryChunkStore` (`ChunkStore.cs`) — unloaded chunks are snapshotted
  (`ChunkSnapshot`) and restored exactly on reload, so relaxed state survives (in-memory now;
  swap for a disk backend later).
- `TerrainStreamer` — a single background worker that owns the model and serves **multiple
  consumers**. Each consumer (`RegisterConsumer`) declares a desired set of coords via
  `SetDesired` and drains ready chunks via `TryDrainResult`. Service levels: *render-ready*
  (chunks + 1-ring halo, seams relaxed, published as deep copies) vs *primal-only* (raw,
  no halo/relax — for gameplay/occupancy). The worker loads the union of all consumers' sets;
  a chunk stays loaded while any consumer wants it. Generation runs in parallel
  (`Parallel.For`); relaxation is incremental (only freshly-generated chunks + neighbours);
  restored/settled chunks are never re-relaxed.
- `ChunkRenderer` — render layer (a render-ready consumer's view): turns published primal
  copies into dual grids + meshes, version-gated. Reads a main-thread *mirror* `TerrainModel`,
  never the worker's authoritative one.
- `ChunkManager` — the MonoBehaviour orchestrator: one render-ready consumer that declares
  `renderRadius` around the camera, installs published copies into the mirror, and meshes.

Debug: `showPrimalGizmos` / `showDualGizmos` on `ChunkManager` and `GridGenerator`.

## Key Constants

- Water density: 1000 kg/m³, viscosity: 1e-3 Pa·s
- Gravity: 9.81 m/s²
- Spatial grid: 5×5×5 cells
- FFT cascade texture: 256×256 (physics readback on cascade 0)
