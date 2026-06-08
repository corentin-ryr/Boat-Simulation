using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Bakery — workshop-scale building wedge + a small square chimney
    // box rising from the cell centroid. Same wedge / wall recipe as
    // BuildingMesh and HouseMesh, sitting at an intermediate roof
    // height. The chimney is the distinguishing silhouette — visible
    // from across the village so the player can find the bakery
    // without inspecting the palette colour.
    //
    // Bakery is the second sub-TileKind of "Building", paired with
    // House. Gameplay: the BakeBread action targets the Baker's
    // assigned Bakery cell; the BuildingRegistry tracks a per-cell
    // Wheat / Bread stockpile.
    public static class BakeryMesh
    {
        const float RoofHeight = 2.6f;
        const float ChimneyHalfWidth = 0.18f;
        const float ChimneyHeight = 0.8f;

        public static Mesh Build(PrimalChunk primal, Func<ChunkCoord, PrimalChunk> lookupPrimal)
        {
            _ = lookupPrimal;

            var verts = new List<Vector3>();
            var tris = new List<int>();

            List<Polygon> dualPolys = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            Vector3[] centroids = primal.PrimalCellCentroids;
            if (dualPolys == null || tp == null || centroids == null) return MeshPrimitives.EmptyMesh();

            int cornerCursor = 0;
            foreach (Polygon poly in dualPolys)
            {
                Vector3[] pverts = poly.GetVerticesPosition();
                int n = pverts.Length;
                if (n < 3) { cornerCursor += n; continue; }

                bool[] isBakery = new bool[n];
                float[] roofY = new float[n];
                bool anyBakery = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = dvpi != null && cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool b = pi >= 0 && pi < tp.Length && tp[pi].Kind == TileKind.Bakery;
                    isBakery[i] = b;
                    if (b) { roofY[i] = centroids[pi].y + RoofHeight; anyBakery = true; }
                }
                cornerCursor += n;
                if (!anyBakery) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isBakery[i]) continue;
                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 ci = pverts[i];
                    Vector3 Mprev = 0.5f * (pverts[iPrev] + ci);
                    Vector3 Mnext = 0.5f * (ci + pverts[iNext]);

                    float ry = roofY[i];
                    Vector3 ciHi = new Vector3(ci.x, ry, ci.z);
                    Vector3 MprevHi = new Vector3(Mprev.x, ry, Mprev.z);
                    Vector3 MnextHi = new Vector3(Mnext.x, ry, Mnext.z);
                    Vector3 Vhi = new Vector3(V.x, ry, V.z);

                    MeshPrimitives.AddTriangle(verts, tris, ciHi, Vhi, MnextHi);
                    MeshPrimitives.AddTriangle(verts, tris, ciHi, MprevHi, Vhi);

                    {
                        bool nbBakery = isBakery[iPrev];
                        float bottomAtV = nbBakery ? roofY[iPrev] : V.y;
                        float bottomAtM = nbBakery ? roofY[iPrev] : Mprev.y;
                        if (ry > bottomAtV || ry > bottomAtM)
                        {
                            Vector3 VlowP = new Vector3(V.x, bottomAtV, V.z);
                            Vector3 MprevLow = new Vector3(Mprev.x, bottomAtM, Mprev.z);
                            MeshPrimitives.EmitWallOutward(verts, tris, VlowP, MprevLow, MprevHi, Vhi, ci);
                        }
                    }
                    {
                        bool nbBakery = isBakery[iNext];
                        float bottomAtV = nbBakery ? roofY[iNext] : V.y;
                        float bottomAtM = nbBakery ? roofY[iNext] : Mnext.y;
                        if (ry > bottomAtV || ry > bottomAtM)
                        {
                            Vector3 MnextLow = new Vector3(Mnext.x, bottomAtM, Mnext.z);
                            Vector3 VlowN = new Vector3(V.x, bottomAtV, V.z);
                            MeshPrimitives.EmitWallOutward(verts, tris, MnextLow, VlowN, Vhi, MnextHi, ci);
                        }
                    }
                }
            }

            // Chimney pass — one small axis-aligned box per Bakery cell,
            // anchored on the cell centroid and rising above the roof.
            // Separate pass (after the wedges) so we don't re-derive
            // roofY per-corner inside the chimney code.
            for (int i = 0; i < tp.Length && i < centroids.Length; i++)
            {
                if (tp[i].Kind != TileKind.Bakery) continue;
                Vector3 c = centroids[i];
                float baseY = c.y + RoofHeight;
                Vector3 min = new Vector3(c.x - ChimneyHalfWidth, baseY,                  c.z - ChimneyHalfWidth);
                Vector3 max = new Vector3(c.x + ChimneyHalfWidth, baseY + ChimneyHeight,  c.z + ChimneyHalfWidth);
                MeshPrimitives.AddBox(verts, tris, min, max);
            }

            return MeshPrimitives.Finalize(verts, tris);
        }
    }
}
