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

## Key Constants

- Water density: 1000 kg/m³, viscosity: 1e-3 Pa·s
- Gravity: 9.81 m/s²
- Spatial grid: 5×5×5 cells
- FFT cascade texture: 256×256 (physics readback on cascade 0)
