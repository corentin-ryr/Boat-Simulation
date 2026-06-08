using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // House — a cottage-scale residential wedge structure. Same geometry
    // recipe as BuildingMesh (raised wedge roof + outward walls on
    // non-same-kind boundaries, terraced merging with neighbours), just
    // a shorter roof height so the silhouette reads as residential next
    // to taller commercial / agricultural Buildings.
    //
    // House is a sub-TileKind of "Building" in player intent: same paint
    // surface (covers the whole cell), same carve behaviour. The
    // distinction is gameplay — an agent's Home is a House cell, and the
    // agent system targets House cells specifically when planning Sleep.
    //
    // Adjacent House cells merge with no internal wall when their roofY
    // matches, same trick BuildingMesh uses. Different-kind adjacency
    // (House ↔ Bakery, House ↔ Building) keeps a wall — the Bakery's
    // taller roof becomes a stepped-back overhang.
    public static class HouseMesh
    {
        const float RoofHeight = 2.0f;

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

                bool[] isHouse = new bool[n];
                float[] roofY = new float[n];
                bool anyHouse = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = dvpi != null && cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool h = pi >= 0 && pi < tp.Length && tp[pi].Kind == TileKind.House;
                    isHouse[i] = h;
                    if (h) { roofY[i] = centroids[pi].y + RoofHeight; anyHouse = true; }
                }
                cornerCursor += n;
                if (!anyHouse) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isHouse[i]) continue;
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

                    // Same wall recipe as BuildingMesh — neighbour-aware
                    // bottom Y so House ↔ House merges flat, House ↔
                    // other steps down to terrain at the boundary.
                    {
                        bool nbHouse = isHouse[iPrev];
                        float bottomAtV = nbHouse ? roofY[iPrev] : V.y;
                        float bottomAtM = nbHouse ? roofY[iPrev] : Mprev.y;
                        if (ry > bottomAtV || ry > bottomAtM)
                        {
                            Vector3 VlowP = new Vector3(V.x, bottomAtV, V.z);
                            Vector3 MprevLow = new Vector3(Mprev.x, bottomAtM, Mprev.z);
                            MeshPrimitives.EmitWallOutward(verts, tris, VlowP, MprevLow, MprevHi, Vhi, ci);
                        }
                    }
                    {
                        bool nbHouse = isHouse[iNext];
                        float bottomAtV = nbHouse ? roofY[iNext] : V.y;
                        float bottomAtM = nbHouse ? roofY[iNext] : Mnext.y;
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
