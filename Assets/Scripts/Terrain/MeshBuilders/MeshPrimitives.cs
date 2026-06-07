using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Small library of mesh-construction helpers shared by every per-kind
    // Build*Mesh class. Kept procedural-only and allocation-modest: each helper
    // appends triangles to the caller's List<Vector3>/List<int> rather than
    // creating intermediate arrays. Finalize() does the one Mesh allocation
    // per chunk per kind.
    //
    // Conventions:
    //   • Winding is counter-clockwise as seen from the visible side, so
    //     RecalculateNormals() produces an outward-facing surface under
    //     Unity's left-handed front-face rule.
    //   • Vertex positions are world-space — the per-chunk parent GameObject
    //     sits at world origin (the dual mesh uses world-Y from primal
    //     vertices) so no transform compensation is needed.
    //   • Empty meshes are returned via EmptyMesh() rather than null so the
    //     mesh assignment / collider cook paths in ChunkSurface never branch
    //     on null and pick up a valid (zero-triangle) Mesh.
    public static class MeshPrimitives
    {
        // One-stop "either empty or populated" closer. Empty path skips the
        // copy entirely; populated path packs the lists into arrays and
        // recalculates normals/bounds.
        public static Mesh Finalize(List<Vector3> vertices, List<int> triangles)
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

        // Stable empty-mesh sentinel — the surface assigns this to MeshFilters
        // and MeshColliders when a kind has zero cells in the chunk, then
        // tears the GameObject down. Marked UInt32 so it can be replaced by
        // a populated mesh of any size without re-cooking the format flag.
        public static Mesh EmptyMesh()
        {
            return new Mesh
            {
                indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
                vertices = System.Array.Empty<Vector3>(),
                triangles = System.Array.Empty<int>(),
            };
        }

        public static void AddTriangle(List<Vector3> verts, List<int> tris,
                                       Vector3 a, Vector3 b, Vector3 c)
        {
            int i = verts.Count;
            verts.Add(a); verts.Add(b); verts.Add(c);
            tris.Add(i); tris.Add(i + 1); tris.Add(i + 2);
        }

        // Convex quad (a, b, c, d) emitted as two tris with shared diagonal a→c.
        // Caller supplies corners in CCW order from the visible side.
        public static void AddQuad(List<Vector3> verts, List<int> tris,
                                   Vector3 a, Vector3 b, Vector3 c, Vector3 d)
        {
            int i = verts.Count;
            verts.Add(a); verts.Add(b); verts.Add(c); verts.Add(d);
            tris.Add(i); tris.Add(i + 1); tris.Add(i + 2);
            tris.Add(i); tris.Add(i + 2); tris.Add(i + 3);
        }

        // World-space axis-aligned box from min to max. Six outward-facing quads.
        // Used by the canopy / tower base primitives.
        public static void AddBox(List<Vector3> v, List<int> t, Vector3 min, Vector3 max)
        {
            Vector3 v0 = new Vector3(min.x, min.y, min.z);
            Vector3 v1 = new Vector3(max.x, min.y, min.z);
            Vector3 v2 = new Vector3(max.x, min.y, max.z);
            Vector3 v3 = new Vector3(min.x, min.y, max.z);
            Vector3 v4 = new Vector3(min.x, max.y, min.z);
            Vector3 v5 = new Vector3(max.x, max.y, min.z);
            Vector3 v6 = new Vector3(max.x, max.y, max.z);
            Vector3 v7 = new Vector3(min.x, max.y, max.z);
            AddQuad(v, t, v4, v5, v6, v7); // +Y top
            AddQuad(v, t, v3, v2, v1, v0); // -Y bottom
            AddQuad(v, t, v0, v1, v5, v4); // -Z
            AddQuad(v, t, v2, v3, v7, v6); // +Z
            AddQuad(v, t, v3, v0, v4, v7); // -X
            AddQuad(v, t, v1, v2, v6, v5); // +X
        }

        // Tapered prism around the vertical axis through (cx, cz): bottom ring
        // at yB with radius rB, top ring at yT with radius rT, sides segments.
        // Watertight (side wall + bottom + top caps). Used for the lighthouse
        // tower stack.
        public static void AddTaperedPrism(List<Vector3> v, List<int> t,
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
            // Side walls. Winding (b0, tB, b1) / (b0, tA, tB) so the
            // triangle normal computed by Cross(p1-p0, p2-p0) points
            // radially outward at the midpoint — matches the
            // wedge convention used elsewhere in MeshBuilders (front face
            // is where the cross product points).
            for (int s = 0; s < sides; s++)
            {
                int s2 = (s + 1) % sides;
                int b0 = baseV + s, b1 = baseV + s2;
                int tA = baseV + sides + s, tB = baseV + sides + s2;
                t.Add(b0); t.Add(tB); t.Add(b1);
                t.Add(b0); t.Add(tA); t.Add(tB);
            }
            // Bottom cap fan — normal points -Y (downward).
            int botC = v.Count;
            v.Add(new Vector3(cx, yB, cz));
            for (int s = 0; s < sides; s++)
            {
                int s2 = (s + 1) % sides;
                t.Add(botC); t.Add(baseV + s); t.Add(baseV + s2);
            }
            // Top cap fan — normal points +Y (upward).
            int topC = v.Count;
            v.Add(new Vector3(cx, yT, cz));
            for (int s = 0; s < sides; s++)
            {
                int s2 = (s + 1) % sides;
                t.Add(topC); t.Add(baseV + sides + s2); t.Add(baseV + sides + s);
            }
        }

        // Vertical wall facing outward (away from `inwardRef`, the wedge's
        // anchor corner). Single-sided — the inside of a building / dock is
        // either the carved-terrain void (no back-face geometry needed) or
        // adjacent same-kind wedge geometry that would z-fight a back face.
        // Picks winding by the sign of the cross product against the outward
        // direction in XZ.
        public static void EmitWallOutward(List<Vector3> verts, List<int> tris,
                                           Vector3 aLow, Vector3 bLow, Vector3 bHigh, Vector3 aHigh,
                                           Vector3 inwardRef)
        {
            Vector3 boundaryMidXZ = 0.5f * (aLow + bLow); boundaryMidXZ.y = 0f;
            Vector3 inwardXZ = inwardRef; inwardXZ.y = 0f;
            Vector3 outwardXZ = boundaryMidXZ - inwardXZ;
            Vector3 edgeDir = bLow - aLow;
            Vector3 upDir = aHigh - aLow;
            Vector3 candidateNormal = Vector3.Cross(edgeDir, upDir);
            candidateNormal.y = 0f;
            if (Vector3.Dot(candidateNormal, outwardXZ) > 0f)
                AddQuad(verts, tris, aLow, bLow, bHigh, aHigh);
            else
                AddQuad(verts, tris, bLow, aLow, aHigh, bHigh);
        }

        // Double-sided vertical wall — used by the road curb where both faces
        // are visible (player walks on either side of the curb).
        public static void EmitWallDoubleSided(List<Vector3> verts, List<int> tris,
                                               Vector3 aLow, Vector3 bLow,
                                               Vector3 bHigh, Vector3 aHigh)
        {
            AddQuad(verts, tris, aLow, bLow, bHigh, aHigh);
            AddQuad(verts, tris, bLow, aLow, aHigh, bHigh);
        }
    }
}
