# Terrain Generation System

## Overview

The terrain is built from a procedural polygon grid using a primal/dual mesh approach. The world is divided into hex-addressed chunks that are streamed in and out as the camera moves. Generation is split into stages that can run asynchronously, with visual quality progressively improving as neighboring chunks become available.

The core pipeline per chunk is:

```
Hex triangles
    → Neighbor computation (edge marking)
    → Random triangle merge (interior only)
    → Subdivision into quads
    → Interior relaxation
    → [Border relaxation — deferred until neighbor is ready]
            → Dual grid computation (read-only cross-chunk lookup at borders)
                    → Mesh build
```

---

## Coordinate System

### `ChunkCoord` — axial hex coordinates `(q, r)`

Chunks are addressed on a hex lattice. The basis vectors in `WorldCenter()` are derived from the inner hex cell geometry so that **border vertices of adjacent chunks always land at exactly the same world positions**, without any negotiation between chunks. This is the invariant the entire border-consistency strategy depends on.

```
WorldCenter(hexRadius R, chunkGridSize N):
    stride = 2 * N
    x = R * 1.5 * stride * q
    z = R * sqrt(3) * stride * (r + q * 0.5)
```

`FromWorldPos` inverts this with standard axial math and cube-coordinate rounding. `HexesInRange(n)` enumerates all chunk coords within `n` hex-rings — used by the streaming manager.

---

## Chunk State Machine

Each chunk progresses through these states. Only `MeshReady` produces visible geometry; earlier states represent work that can be done off the main thread.

```
Unloaded
    → PrimalReady       (hex gen + merge + subdivide, async)
    → InteriorRelaxed   (relax non-edge vertices, async)
    → BorderRelaxed     (joint pass with each loaded neighbor, async)
    → DualReady         (dual computation, async)
    → MeshReady         (mesh upload, main thread)
```

A chunk can reach `DualReady` without all neighbors being `InteriorRelaxed` — border vertices will simply be at their unrelaxed hex grid positions for the missing sides. When a new neighbor finishes `InteriorRelaxed`, the border pass fires for both chunks and both advance.

---

## Primal Grid Pipeline

All steps here are per-chunk with no cross-chunk data. They can run entirely on a background thread.

### Step 1 — Hex shape generation (`GenHexShape`)

A hex-shaped region of `chunkGridSize` rings is laid out in axial coordinates. Each hex cell is split into 6 triangles (center + 6 corner vertices). Vertices are stored in a `VertexCollection` keyed by rounded world position to deduplicate shared corners.

No ghost ring. The grid is exactly `chunkGridSize` rings, no more.

### Step 2 — Neighbor computation and edge marking (`ComputeNeighbors`)

For each triangle, shared edges are found by checking vertex pairs. Any triangle edge with no neighbor gets its two vertices flagged as `IsEdge = true`.

`IsEdge` vertices are the chunk's border. They sit at fixed world positions determined only by `hexRadius`, `chunkGridSize`, and the chunk's `WorldCenter()` — no random element, no dependence on other chunks.

### Step 3 — Random triangle merge (`RandomTriangleMerge`)

Adjacent interior triangle pairs are randomly merged into quads. The constraint: **a triangle is never merged if any of its vertices is `IsEdge`**. This keeps border triangle topology deterministic, which is what guarantees that adjacent chunks' border vertices are at identical world positions.

Implementation note: the current code uses `List<Polygon>.Remove()` which is O(P). Replace with a swap-and-pop pattern (swap the target with the last element, then remove the last) to make the loop O(P) total instead of O(P²).

### Step 4 — Polygon subdivision (`SubdividePolygon`)

Each triangle and quad is subdivided by inserting a center vertex and edge midpoints, producing a grid of quads. After this step all polygons have exactly 4 vertices.

### Step 5 — Interior relaxation (`GridRelaxation`)

A Lloyd-style relaxation pass runs for `nbIterRelaxation` iterations over **non-edge vertices only**. Each quad accumulates a rotational goal (aligning opposite vertices 90° apart) and vertices move 10% toward it per iteration.

`IsEdge` vertices do not move. They remain at their deterministic hex grid positions, acting as fixed anchors for the relaxation. This ensures the border vertex positions are unchanged after this step.

---

## Border Relaxation

### Motivation

After interior relaxation, each chunk's border vertices are still at their initial hex grid positions. The quads adjacent to the border are therefore not as regular as the interior. When a neighboring chunk is available, a joint relaxation pass can move border vertices using the full quad neighborhood from both sides, making the boundary region as regular as the interior.

### Triggering condition

The border relaxation between chunk A and chunk B fires when **both** chunks have reached `InteriorRelaxed`. It is a one-time pass per pair — once border vertices have been moved, the result is committed and no further border relaxation is run for that edge.

### Algorithm

For a given shared edge between chunk A and chunk B:

1. **Find matching border vertex pairs.** Border vertices in both chunks are at identical world positions (guaranteed by the coordinate system). Match them by rounding their positions to 3 decimal places — the same rounding used by `VertexCollection.Key()`.

2. **Gather the full quad neighborhood.** For each matched pair (V_A, V_B), collect:
   - All quads in chunk A that contain V_A
   - All quads in chunk B that contain V_B
   This is the complete local neighborhood that would surround V if the chunks were stitched.

3. **Compute the joint relaxation goal.** Run the standard relaxation formula over the combined quad set — the same formula used in `GridRelaxation`, treating all quads as if they belong to one mesh.

4. **Apply the same movement to both copies.** The resulting delta is applied identically to V_A and to V_B. After this step both copies are at the same new world position.

No data structures are merged. Each chunk remains self-contained. The pass is purely a read from one chunk's polygon list while writing to both chunks' vertex positions.

### Corner vertices (three-chunk intersections)

Vertices at the corner of three chunks need special handling: the joint pass for the A-B edge should not move the corner vertex yet if chunk C (the third neighbor) is not yet `InteriorRelaxed`. Defer the corner vertex movement until all three neighbors are ready. Corner vertices can be identified as those that appear in more than one shared edge.

---

## Dual Grid Computation

The dual grid is computed after a chunk reaches `BorderRelaxed` (or `InteriorRelaxed` if some neighbors are missing). It produces the final polygon shapes that become the terrain mesh.

### Standard dual (interior primal vertices)

For each non-edge primal vertex V, create a dual polygon whose vertices are the centers of all quads surrounding V, ordered by signed angle. This is the existing `GenerateDualGrid` logic.

### Border dual (edge primal vertices)

For each edge primal vertex V, the full quad neighborhood spans two chunks. The border dual polygon for V is built by:

1. Collecting quad centers from V's own chunk (quads in the local `VertexCollection`).
2. **Read-only cross-chunk lookup**: find the matching border vertex in each loaded neighbor and collect their quad centers.
3. Sort all centers by signed angle around V's world position and build the dual polygon.

This is not stitching — neither chunk's `VertexCollection` is modified. The neighbor's data is only read to complete the polygon.

If a neighbor is not yet loaded, its quads are simply absent. The border dual polygon is built from whatever is available and will be rebuilt when the neighbor loads and the border relaxation fires.

### Output

The dual computation produces a `List<Polygon>` for the chunk. Only polygons whose center maps back to this chunk's `ChunkCoord` are kept (same filter as before, but now with correct cross-chunk vertices so the filter behaves consistently).

---

## Mesh Building

Unchanged from the current implementation: fan-triangulate each dual polygon and upload as a Unity `Mesh` with `IndexFormat.UInt32`. Mesh upload must happen on the main thread; everything before it can be async.

When border relaxation fires for an already-meshed chunk, the dual and mesh for the affected boundary strip must be recomputed. The simplest approach is to rebuild the entire chunk mesh — dual computation is cheap (O(P)) and mesh upload is fast for chunk-sized meshes.

---

## Streaming — `TerrainStreamer` + `ChunkManager`

### Layers

- **`TerrainModel`** — authoritative, Unity-free data: loaded `PrimalChunk`s, generation, joint
  border relaxation, and the chunk store. Safe to run off the main thread.
- **`TerrainStreamer`** — a single dedicated background thread that owns the model and serves
  multiple **consumers**. A consumer (`RegisterConsumer(renderReady)`) declares a desired set of
  coords with `SetDesired(consumer, coords)` and drains ready chunks with
  `TryDrainResult(consumer, …)`. The worker loads the **union** of all consumers' effective sets;
  a chunk is unloaded only when no consumer wants it.
- **`ChunkRenderer`** — turns published primal copies into dual grids + meshes (version-gated),
  reading a main-thread **mirror** `TerrainModel`, never the worker's.
- **`ChunkManager`** — MonoBehaviour orchestrator: one render-ready consumer that declares
  `renderRadius` around the camera each time it crosses a chunk boundary.

### Service levels

- **render-ready** — the consumer's set is expanded by a **1-ring halo** (so seams can weld and
  the renderer's cross-chunk dual lookup has neighbours), the seams are relaxed, and chunks are
  published as independent **deep copies** (`PrimalChunk.DeepCopy`) for the main thread to install.
- **primal-only** — chunks are loaded as generated (borders pinned), with **no halo and no
  relaxation**; isolated/non-contiguous requests are valid. For gameplay/occupancy/AI.

### Worker pass

Per pass the worker: evicts chunks no consumer needs (persisting via the store), restores stored
chunks serially, **generates fresh chunks in parallel** (`Parallel.For` — each is independent),
runs **incremental** border relaxation (only freshly-generated render-region chunks + their
neighbours; restored/settled chunks are never re-relaxed), then publishes per consumer the
version-changed chunks in its effective set. Camera moves coalesce (newest desired set wins); a
pass abandons early if a newer request arrives mid-flight.

### Persistence

Unloaded chunks are snapshotted (`ChunkSnapshot` → `IChunkStore`) and restored **exactly** on
reload, so a chunk's relaxed/welded state survives instead of regenerating fresh. `MemoryChunkStore`
is per-session; a disk backend can slot in behind `IChunkStore`.

### Chunk seed

Each chunk's random seed is derived deterministically:

```csharp
int chunkSeed = globalSeed ^ (coord.q * 73856093) ^ (coord.r * 19349663);
```

This ensures the same chunk always generates the same primal topology regardless of load order.

---

## What this design does NOT do

- **No ghost ring.** Primal grids are generated at their natural size (`chunkGridSize` rings). The ghost ring was an approximation that produced geometric cracks; it is replaced by joint cross-chunk border relaxation and cross-chunk dual completion (both implemented, keyed on `LatticeKey`).
- **No stitching.** `VertexCollection`s are never merged. Each chunk owns its own primal data permanently; seams stay welded because shared border vertices are moved jointly during relaxation.
- **Border triangles DO merge.** The merge pass is allowed to reach the seams (no border guard), so randomness isn't cut off into a regular band. Boundary edges have no neighbour and never merge, so the seam lattice points stay deterministic and coincide across chunks regardless of how interior merging resolves on each side.
- **No LOD.** All chunks are generated at the same resolution. A LOD system would be a separate concern layered on top of this pipeline.
