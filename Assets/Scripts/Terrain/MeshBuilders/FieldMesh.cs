using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Flat tinted wedge sheet covering each Field cell, lifted by a hair to
    // win the z-fight at the carved-wedge boundary. The terrain mesh side
    // skips the wedge of every Field corner (Field is in the carve list),
    // so this overlay IS the visible cell surface — not a decal on top.
    //
    // No walls, no curbs — the cell sits flush with neighbour terrain. The
    // Field material's base colour carries the visual identity (green crop
    // rows). For a future Plaza / Sand kind, copy this file and tweak `Lift`
    // + material wiring; the rest of the logic is identical.
    public static class FieldMesh
    {
        const float Lift = 0.03f;

        public static Mesh Build(PrimalChunk primal, Func<ChunkCoord, PrimalChunk> lookupPrimal)
        {
            _ = lookupPrimal;

            var vertices = new List<Vector3>();
            var triangles = new List<int>();

            List<Polygon> polygons = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            if (polygons == null || dvpi == null || tp == null) return MeshPrimitives.EmptyMesh();

            Vector3 liftV = new Vector3(0f, Lift, 0f);

            int cornerCursor = 0;
            foreach (Polygon poly in polygons)
            {
                Vector3[] verts = poly.GetVerticesPosition();
                int n = verts.Length;
                if (n < 3) { cornerCursor += n; continue; }

                bool[] isField = new bool[n];
                bool any = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool k = pi >= 0 && pi < tp.Length && tp[pi].Kind == TileKind.Field;
                    isField[i] = k;
                    if (k) any = true;
                }
                cornerCursor += n;
                if (!any) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += verts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isField[i]) continue;
                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 Mprev = 0.5f * (verts[iPrev] + verts[i]);
                    Vector3 Mnext = 0.5f * (verts[i] + verts[iNext]);

                    // (c_i, V, M_next) and (c_i, M_prev, V) — same winding
                    // as the terrain mesh's wedge geometry.
                    int b = vertices.Count;
                    vertices.Add(verts[i] + liftV);
                    vertices.Add(Mnext + liftV);
                    vertices.Add(V + liftV);
                    vertices.Add(Mprev + liftV);
                    triangles.Add(b + 0); triangles.Add(b + 2); triangles.Add(b + 1);
                    triangles.Add(b + 0); triangles.Add(b + 3); triangles.Add(b + 2);
                }
            }

            return MeshPrimitives.Finalize(vertices, triangles);
        }
    }
}
