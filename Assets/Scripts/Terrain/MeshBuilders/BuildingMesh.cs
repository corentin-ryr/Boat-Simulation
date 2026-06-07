using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Raised-roof building overlay: each Building cell's wedge is lifted to
    // a per-cell roof height (`primalCentroid.y + BuildingHeight`), and
    // vertical walls descend from the roof to the terrain (or to the
    // lower-roofed neighbour's roof, for adjacent buildings with stepped
    // heights). Adjacent Building cells with the same roofY merge into one
    // flat roof with no internal wall.
    //
    // Terrain BuildMesh has already carved away the Building wedges, so
    // this overlay is filling the matching hole (and adding walls /
    // overhangs).
    public static class BuildingMesh
    {
        const float BuildingHeight = 3.0f;

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

                bool[] isBuilding = new bool[n];
                float[] roofY = new float[n];
                bool anyBuilding = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = dvpi != null && cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool b = pi >= 0 && pi < tp.Length && tp[pi].Kind == TileKind.Building;
                    isBuilding[i] = b;
                    if (b) { roofY[i] = centroids[pi].y + BuildingHeight; anyBuilding = true; }
                }
                cornerCursor += n;
                if (!anyBuilding) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isBuilding[i]) continue;
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

                    // V→Mprev wall. Neighbour wedge iPrev controls the
                    // bottom Y: Building → its flat roofY (constant);
                    // non-Building → per-endpoint terrain Y so the wall
                    // tracks the slope.
                    {
                        bool nbBuilding = isBuilding[iPrev];
                        float bottomAtV = nbBuilding ? roofY[iPrev] : V.y;
                        float bottomAtM = nbBuilding ? roofY[iPrev] : Mprev.y;
                        if (ry > bottomAtV || ry > bottomAtM)
                        {
                            Vector3 VlowP = new Vector3(V.x, bottomAtV, V.z);
                            Vector3 MprevLow = new Vector3(Mprev.x, bottomAtM, Mprev.z);
                            MeshPrimitives.EmitWallOutward(verts, tris, VlowP, MprevLow, MprevHi, Vhi, ci);
                        }
                    }
                    {
                        bool nbBuilding = isBuilding[iNext];
                        float bottomAtV = nbBuilding ? roofY[iNext] : V.y;
                        float bottomAtM = nbBuilding ? roofY[iNext] : Mnext.y;
                        if (ry > bottomAtV || ry > bottomAtM)
                        {
                            Vector3 MnextLow = new Vector3(Mnext.x, bottomAtM, Mnext.z);
                            Vector3 VlowN = new Vector3(V.x, bottomAtV, V.z);
                            MeshPrimitives.EmitWallOutward(verts, tris, MnextLow, VlowN, Vhi, MnextHi, ci);
                        }
                    }
                }
            }

            return MeshPrimitives.Finalize(verts, tris);
        }
    }
}
