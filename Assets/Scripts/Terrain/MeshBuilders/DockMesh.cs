using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Harbour platform — a raised stone wedge covering every Dock OR
    // Lighthouse cell. Renamed conceptually from "Dock" because the
    // Lighthouse shares the exact same base; only the tower on top
    // distinguishes one from the other. Owned by the Dock overlay so
    // there is a single material for the whole harbour.
    //
    // Geometry per platform corner of a dual polygon:
    //   • Flat top wedge (c_i, V, M_next, M_prev) at fixed world Y =
    //     PlatformTopY. The deck never tilts with the seabed; deeper
    //     water just produces a taller wall.
    //   • Vertical outward-facing wall on every wedge boundary where the
    //     neighbour wedge isn't also a platform (Dock or Lighthouse).
    //     Single-sided. Wall base Y = min(V.y, neighbour.y) - BuriedDepth
    //     so it sinks into seabed / slope and never floats.
    //
    // Adjacent platform cells merge with no internal wall — Dock↔Dock,
    // Dock↔Lighthouse, Lighthouse↔Lighthouse all read as one continuous
    // embankment.
    public static class DockMesh
    {
        // Sea level convention used across the boat-sim.
        const float SeaLevelY = 0f;
        // Absolute world Y of the harbour deck. Lighthouse tower base
        // must match this value (LighthouseMesh.PlatformTopY).
        const float PlatformTopY = SeaLevelY + 2f;
        // How far the outer wall buries into the terrain below the wedge
        // boundary, so sloped ground never reveals a gap underneath.
        const float BuriedDepth = 0.1f;

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

                bool[] isPlatform = new bool[n];
                bool anyPlatform = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = dvpi != null && cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    TileKind k = (pi >= 0 && pi < tp.Length) ? tp[pi].Kind : TileKind.Default;
                    bool p = (k == TileKind.Dock || k == TileKind.Lighthouse);
                    isPlatform[i] = p;
                    if (p) anyPlatform = true;
                }
                cornerCursor += n;
                if (!anyPlatform) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isPlatform[i]) continue;
                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 ci = pverts[i];
                    Vector3 Mprev = 0.5f * (pverts[iPrev] + ci);
                    Vector3 Mnext = 0.5f * (ci + pverts[iNext]);

                    // Deck top sits at a fixed world Y, independent of
                    // the cell's centroid / seabed depth. This is the
                    // whole point of the "uniform deck" rule.
                    float ry = PlatformTopY;
                    Vector3 ciHi = new Vector3(ci.x, ry, ci.z);
                    Vector3 MprevHi = new Vector3(Mprev.x, ry, Mprev.z);
                    Vector3 MnextHi = new Vector3(Mnext.x, ry, Mnext.z);
                    Vector3 Vhi = new Vector3(V.x, ry, V.z);

                    // Flat top wedge — (c_i, V, M_next) / (c_i, M_prev,
                    // V) winding matches BuildingMesh so normals point up.
                    MeshPrimitives.AddTriangle(verts, tris, ciHi, Vhi, MnextHi);
                    MeshPrimitives.AddTriangle(verts, tris, ciHi, MprevHi, Vhi);

                    // V→Mprev wall: emit when the neighbour wedge isn't
                    // also a platform. Buried base so sloped terrain
                    // doesn't reveal a floating gap underneath.
                    if (!isPlatform[iPrev])
                    {
                        float baseY = Mathf.Min(V.y, Mprev.y) - BuriedDepth;
                        Vector3 VlowP = new Vector3(V.x, baseY, V.z);
                        Vector3 MprevLow = new Vector3(Mprev.x, baseY, Mprev.z);
                        MeshPrimitives.EmitWallOutward(verts, tris, VlowP, MprevLow, MprevHi, Vhi, ci);
                    }
                    if (!isPlatform[iNext])
                    {
                        float baseY = Mathf.Min(V.y, Mnext.y) - BuriedDepth;
                        Vector3 MnextLow = new Vector3(Mnext.x, baseY, Mnext.z);
                        Vector3 VlowN = new Vector3(V.x, baseY, V.z);
                        MeshPrimitives.EmitWallOutward(verts, tris, MnextLow, VlowN, Vhi, MnextHi, ci);
                    }
                }
            }

            return MeshPrimitives.Finalize(verts, tris);
        }
    }
}
