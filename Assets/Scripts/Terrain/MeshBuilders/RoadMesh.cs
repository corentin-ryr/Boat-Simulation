using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Townscaper-style road overlay: each Road cell's wedge is sunk by
    // RoadSinkDepth (so the surface sits below the terrain that's been carved
    // away for it), and any wedge boundary that touches a non-Road wedge of
    // the same dual polygon gets a vertical curb running from the sunk Y up
    // to terrain Y. Adjacent Road wedges share their midpoint at the same
    // sunken height, so the road surface is seamless across cell joins.
    //
    // Cross-chunk note: dvpi marks foreign corners as -1, so a road that
    // runs across a chunk seam is currently undetected at the foreign end of
    // the shared dual polygon. The cross-chunk Version cascade keeps both
    // sides rebuilt in sync; resolving foreign tile kinds per dual corner
    // would need a data-layer extension and is deferred.
    public static class RoadMesh
    {
        const float RoadSinkDepth = 0.06f;

        // lookupPrimal isn't read today but is part of the standard
        // mesh-builder signature so future seam-aware logic can land without
        // breaking the dispatcher in ChunkSurface.
        public static Mesh Build(PrimalChunk primal, Func<ChunkCoord, PrimalChunk> lookupPrimal)
        {
            _ = lookupPrimal;

            var verts = new List<Vector3>();
            var tris = new List<int>();

            List<Polygon> dualPolys = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            if (dualPolys == null || tp == null) return MeshPrimitives.EmptyMesh();

            Vector3 sink = new Vector3(0, RoadSinkDepth, 0);

            int cornerCursor = 0;
            foreach (Polygon poly in dualPolys)
            {
                Vector3[] pverts = poly.GetVerticesPosition();
                int n = pverts.Length;
                if (n < 3) { cornerCursor += n; continue; }

                bool[] isRoad = new bool[n];
                bool anyRoad = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = dvpi != null && cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool r = pi >= 0 && pi < tp.Length && tp[pi].Kind == TileKind.Road;
                    isRoad[i] = r;
                    if (r) anyRoad = true;
                }
                cornerCursor += n;
                if (!anyRoad) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;
                Vector3 Vlow = V - sink;

                for (int i = 0; i < n; i++)
                {
                    if (!isRoad[i]) continue;
                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 ci = pverts[i];
                    Vector3 Mprev = 0.5f * (pverts[iPrev] + ci);
                    Vector3 Mnext = 0.5f * (ci + pverts[iNext]);

                    Vector3 ciLow = ci - sink;
                    Vector3 MprevLow = Mprev - sink;
                    Vector3 MnextLow = Mnext - sink;

                    // Sunk road surface — winding matches the terrain wedge
                    // it replaces so the normal points up.
                    MeshPrimitives.AddTriangle(verts, tris, ciLow, Vlow, MnextLow);
                    MeshPrimitives.AddTriangle(verts, tris, ciLow, MprevLow, Vlow);

                    // Curbs: V→M boundaries inside the polygon where the
                    // neighbour wedge isn't road. Double-sided so visible
                    // from inside the road and from the grass approach.
                    if (!isRoad[iPrev])
                        MeshPrimitives.EmitWallDoubleSided(verts, tris, Vlow, MprevLow, Mprev, V);
                    if (!isRoad[iNext])
                        MeshPrimitives.EmitWallDoubleSided(verts, tris, MnextLow, Vlow, V, Mnext);
                }
            }

            return MeshPrimitives.Finalize(verts, tris);
        }
    }
}
