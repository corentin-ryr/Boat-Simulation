using System;
using UnityEngine;

namespace TerrainGrid
{
    // A presenter owns how a single TileKind is *visualised* in the scene. It is a
    // ScriptableObject (and so authored, not coded, into the project) so the same
    // C# class can be parameterised through several .asset variants — e.g. one
    // WedgeAreaPresenter asset configured for Field, another for a future Plaza
    // tile, both sharing the same code path with different sink / lift / colour
    // parameters.
    //
    // Lifecycle is split between this abstract class and the concrete subclass:
    //   - The base class' Rebuild() handles the bookkeeping that's identical
    //     across every presenter: lazy-spawn the child GameObject, version-gate
    //     the rebuild, swap meshes safely, and tear the GO down when the chunk
    //     stops having any cells of this kind. Subclasses never see this.
    //   - The subclass overrides BuildMesh() and returns a freshly constructed
    //     UnityEngine.Mesh for the chunk's current cell set. It must NOT mutate
    //     `handle.*` — the base class owns that.
    //
    // Cross-chunk lookups (the road overlay's curb decisions, the building's
    // wall step-cap heights) come in through `PresenterContext.lookupPrimal`,
    // so a presenter can read neighbour chunks' TileProperties without coupling
    // to ChunkManager directly.
    public abstract class TilePresenter : ScriptableObject
    {
        // True if the surface should attach a MeshCollider tracking the overlay
        // mesh. Defaults off — most overlays are decorative (Field is flat on
        // walkable terrain, Market/Lighthouse/Dock are point structures we let
        // the player walk through for v1). Road and Building override to true
        // because their carved-terrain holes leave no ground / wall to block
        // the player without an overlay collider.
        public virtual bool RequiresCollider => false;

        // Drives the once-per-Apply rebuild cycle. Called by ChunkSurface from
        // EnsurePresence whenever a non-flat chunk wants rendering AND has cells
        // of this presenter's kind. Returns whether the chunk actually has any
        // cell of `ctx.kind` (so the caller can detach the collider if it just
        // emptied out).
        public bool Rebuild(in PresenterContext ctx)
        {
            ChunkOverlayHandle handle = ctx.handle;
            if (handle == null) return false;

            bool anyOfKind = HasAnyKind(ctx.primal, ctx.kind);

            if (!anyOfKind)
            {
                // Chunk used to have ≥1 cell of this kind, doesn't now → free the
                // child. Cheap when nothing's there already (Go == null short-circuit).
                if (handle.Go != null) DestroyOverlay(handle);
                handle.BuiltFromVersion = ctx.primal.Version;
                return false;
            }

            // Same Version since the last build → existing Mesh is still correct.
            if (handle.BuiltFromVersion == ctx.primal.Version && handle.Mesh != null) return true;

            if (handle.Go == null)
            {
                handle.Go = new GameObject(ctx.label);
                handle.Go.transform.SetParent(ctx.parent, false);
                handle.Mf = handle.Go.AddComponent<MeshFilter>();
                handle.Mr = handle.Go.AddComponent<MeshRenderer>();
                handle.Mr.sharedMaterial = ctx.material;
            }

            // Free the old Mesh asset before replacing it — sharedMesh = newMesh
            // leaves the old asset orphaned otherwise.
            if (handle.Mesh != null) UnityEngine.Object.Destroy(handle.Mesh);
            handle.Mesh = BuildMesh(in ctx);
            handle.Mesh.name = $"{ctx.label} {ctx.coord}";
            handle.Mf.sharedMesh = handle.Mesh;
            handle.BuiltFromVersion = ctx.primal.Version;
            return true;
        }

        // Subclass entry point. Builds the chunk-combined overlay mesh from the
        // current cell set. Implementations should:
        //   - Walk `ctx.primal.TileProperties` (or `ctx.primal.Dual` for area
        //     presenters) to find cells of `ctx.kind`.
        //   - Emit geometry for each. World-space positions; the overlay GO is
        //     parented at world origin (same as the terrain mesh) so vertex
        //     positions go in as-is.
        //   - Return an empty Mesh (`MeshUtils.EmptyMesh()`) if there's nothing
        //     to draw rather than throwing.
        protected abstract Mesh BuildMesh(in PresenterContext ctx);

        // Tear down the overlay GO + Mesh. Called by Rebuild() when the chunk
        // drops to zero cells of this kind, and by ChunkSurface directly when a
        // chunk leaves the wanted set entirely (eviction path).
        public virtual void DestroyOverlay(ChunkOverlayHandle handle)
        {
            if (handle == null) return;
            if (handle.Mc != null) { UnityEngine.Object.Destroy(handle.Mc); handle.Mc = null; }
            if (handle.Mesh != null) { UnityEngine.Object.Destroy(handle.Mesh); handle.Mesh = null; }
            if (handle.Go != null) { UnityEngine.Object.Destroy(handle.Go); handle.Go = null; }
            handle.Mf = null;
            handle.Mr = null;
            handle.BuiltFromVersion = -1;
        }

        // Common predicate factored here so every concrete presenter doesn't
        // need to re-implement it. Treats a null TileProperties array as
        // "this chunk hasn't been finalised yet" → no cells of any kind.
        protected static bool HasAnyKind(PrimalChunk primal, TileKind kind)
        {
            TileProperty[] tp = primal.TileProperties;
            if (tp == null) return false;
            for (int i = 0; i < tp.Length; i++) if (tp[i].Kind == kind) return true;
            return false;
        }
    }

    // Per-chunk per-kind overlay state. One instance per (Chunk, TileKind) pair
    // that's ever held an overlay, stored in Presence.Overlays. Owned by
    // ChunkSurface; presenters mutate the fields through their Rebuild() /
    // DestroyOverlay() calls, but never construct or store a reference.
    public class ChunkOverlayHandle
    {
        public GameObject Go;
        public MeshFilter Mf;
        public MeshRenderer Mr;
        public Mesh Mesh;
        public MeshCollider Mc;          // populated by ChunkSurface when the presenter has RequiresCollider=true AND wantCollider
        public int BuiltFromVersion = -1;
    }

    // Argument bundle passed through TilePresenter.Rebuild()/BuildMesh(). Struct
    // + `in` so the call site reads naturally without per-Apply allocation. The
    // `lookupPrimal` delegate is the presenter's hook for cross-chunk reads —
    // foreign neighbour cells live in OTHER chunks, which the presenter resolves
    // by axial offset + lookup.
    public readonly struct PresenterContext
    {
        public readonly ChunkCoord coord;
        public readonly PrimalChunk primal;
        public readonly TileKind kind;
        public readonly Transform parent;
        public readonly Material material;
        public readonly ChunkOverlayHandle handle;
        public readonly Func<ChunkCoord, PrimalChunk> lookupPrimal;
        public readonly string label;

        public PresenterContext(ChunkCoord coord, PrimalChunk primal, TileKind kind,
                                Transform parent, Material material,
                                ChunkOverlayHandle handle,
                                Func<ChunkCoord, PrimalChunk> lookupPrimal, string label)
        {
            this.coord = coord;
            this.primal = primal;
            this.kind = kind;
            this.parent = parent;
            this.material = material;
            this.handle = handle;
            this.lookupPrimal = lookupPrimal;
            this.label = label;
        }
    }

    // Small helpers presenters use to assemble Meshes. Public so any new
    // presenter living outside the Terrain namespace can pull them in.
    public static class MeshUtils
    {
        public static Mesh EmptyMesh()
        {
            return new Mesh
            {
                indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
                vertices = Array.Empty<Vector3>(),
                triangles = Array.Empty<int>(),
            };
        }

        public static Mesh Finalize(System.Collections.Generic.List<Vector3> vertices,
                                    System.Collections.Generic.List<int> triangles)
        {
            if (vertices.Count == 0) return EmptyMesh();
            Mesh m = new Mesh
            {
                indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
            };
            m.SetVertices(vertices);
            m.SetTriangles(triangles, 0);
            m.RecalculateNormals();
            m.RecalculateBounds();
            return m;
        }

        public static void AddTriangle(System.Collections.Generic.List<Vector3> verts,
                                       System.Collections.Generic.List<int> tris,
                                       Vector3 a, Vector3 b, Vector3 c)
        {
            int i = verts.Count;
            verts.Add(a); verts.Add(b); verts.Add(c);
            tris.Add(i); tris.Add(i + 1); tris.Add(i + 2);
        }

        public static void AddQuad(System.Collections.Generic.List<Vector3> verts,
                                   System.Collections.Generic.List<int> tris,
                                   Vector3 a, Vector3 b, Vector3 c, Vector3 d)
        {
            int i = verts.Count;
            verts.Add(a); verts.Add(b); verts.Add(c); verts.Add(d);
            tris.Add(i); tris.Add(i + 1); tris.Add(i + 2);
            tris.Add(i); tris.Add(i + 2); tris.Add(i + 3);
        }

        // World-space axis-aligned box. Six outward-facing quads.
        public static void AddBox(System.Collections.Generic.List<Vector3> v,
                                  System.Collections.Generic.List<int> t,
                                  Vector3 min, Vector3 max)
        {
            Vector3 v0 = new Vector3(min.x, min.y, min.z);
            Vector3 v1 = new Vector3(max.x, min.y, min.z);
            Vector3 v2 = new Vector3(max.x, min.y, max.z);
            Vector3 v3 = new Vector3(min.x, min.y, max.z);
            Vector3 v4 = new Vector3(min.x, max.y, min.z);
            Vector3 v5 = new Vector3(max.x, max.y, min.z);
            Vector3 v6 = new Vector3(max.x, max.y, max.z);
            Vector3 v7 = new Vector3(min.x, max.y, max.z);
            AddQuad(v, t, v4, v5, v6, v7);
            AddQuad(v, t, v3, v2, v1, v0);
            AddQuad(v, t, v0, v1, v5, v4);
            AddQuad(v, t, v2, v3, v7, v6);
            AddQuad(v, t, v3, v0, v4, v7);
            AddQuad(v, t, v1, v2, v6, v5);
        }

        // Tapered prism around the vertical axis through (cx, cz). Used for the
        // lighthouse's stacked sections.
        public static void AddTaperedPrism(System.Collections.Generic.List<Vector3> v,
                                           System.Collections.Generic.List<int> t,
                                           float cx, float cz,
                                           float yB, float rB, float yT, float rT, int sides)
        {
            int baseV = v.Count;
            for (int s = 0; s < sides; s++)
            {
                float a = (s / (float)sides) * Mathf.PI * 2f;
                v.Add(new Vector3(cx + rB * Mathf.Cos(a), yB, cz + rB * Mathf.Sin(a)));
            }
            for (int s = 0; s < sides; s++)
            {
                float a = (s / (float)sides) * Mathf.PI * 2f;
                v.Add(new Vector3(cx + rT * Mathf.Cos(a), yT, cz + rT * Mathf.Sin(a)));
            }
            for (int s = 0; s < sides; s++)
            {
                int s2 = (s + 1) % sides;
                int b0 = baseV + s, b1 = baseV + s2;
                int t0 = baseV + sides + s, t1 = baseV + sides + s2;
                t.Add(b0); t.Add(b1); t.Add(t1);
                t.Add(b0); t.Add(t1); t.Add(t0);
            }
            int botC = v.Count;
            v.Add(new Vector3(cx, yB, cz));
            for (int s = 0; s < sides; s++)
            {
                int s2 = (s + 1) % sides;
                t.Add(botC); t.Add(baseV + s2); t.Add(baseV + s);
            }
            int topC = v.Count;
            v.Add(new Vector3(cx, yT, cz));
            for (int s = 0; s < sides; s++)
            {
                int s2 = (s + 1) % sides;
                t.Add(topC); t.Add(baseV + sides + s); t.Add(baseV + sides + s2);
            }
        }
    }
}
