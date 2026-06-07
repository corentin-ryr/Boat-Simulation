using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Market as a wedge-shaped plaza + one centroid-anchored canopy slab
    // per cell. The plaza fills the cell footprint (same wedge fan as Field
    // but raised), making it clear the whole cell is the market. The canopy
    // sits above the centroid to give the silhouette a recognisable
    // "open-air stall" shape from the RTS view.
    //
    // Two passes in one mesh:
    //   1. Plaza wedge fan — every Market corner of every dual polygon.
    //   2. Per-Market-cell canopy — one tilted slab per cell at its centroid,
    //      anchored at canopy bottom = plaza top + CanopyClearance.
    //
    // No carving (Market isn't in IsCarvingKind), no collider for v1 — the
    // player walks under the canopy and onto the plaza freely.
    public static class MarketMesh
    {
        // Plaza raise.
        const float PlazaLift = 0.2f;
        // Canopy geometry (footprint at centroid, ported from old MarketPresenter).
        const float CanopyHalf = 0.30f;
        const float CanopyOverhang = 0.10f;
        const float CanopyClearance = 0.4f;   // gap from plaza top to canopy back edge
        const float CanopySlope = 0.45f;      // front-edge Y = back-edge Y - CanopySlope

        public static Mesh Build(PrimalChunk primal, Func<ChunkCoord, PrimalChunk> lookupPrimal)
        {
            _ = lookupPrimal;

            var verts = new List<Vector3>();
            var tris = new List<int>();

            List<Polygon> polygons = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            Vector3[] centroids = primal.PrimalCellCentroids;
            if (polygons == null || dvpi == null || tp == null || centroids == null)
                return MeshPrimitives.EmptyMesh();

            Vector3 lift = new Vector3(0f, PlazaLift, 0f);

            // --- Pass 1: plaza wedges (per dual polygon corner of Market kind) ---
            int cornerCursor = 0;
            foreach (Polygon poly in polygons)
            {
                Vector3[] pverts = poly.GetVerticesPosition();
                int n = pverts.Length;
                if (n < 3) { cornerCursor += n; continue; }

                bool[] isMarket = new bool[n];
                bool any = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool m = pi >= 0 && pi < tp.Length && tp[pi].Kind == TileKind.Market;
                    isMarket[i] = m;
                    if (m) any = true;
                }
                cornerCursor += n;
                if (!any) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += pverts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isMarket[i]) continue;
                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 Mprev = 0.5f * (pverts[iPrev] + pverts[i]);
                    Vector3 Mnext = 0.5f * (pverts[i] + pverts[iNext]);

                    int b = verts.Count;
                    verts.Add(pverts[i] + lift);
                    verts.Add(Mnext + lift);
                    verts.Add(V + lift);
                    verts.Add(Mprev + lift);
                    tris.Add(b + 0); tris.Add(b + 2); tris.Add(b + 1);
                    tris.Add(b + 0); tris.Add(b + 3); tris.Add(b + 2);
                }
            }

            // --- Pass 2: canopy slabs (one per Market cell at its centroid) ---
            for (int i = 0; i < tp.Length && i < centroids.Length; i++)
            {
                if (tp[i].Kind != TileKind.Market) continue;
                Vector3 c = centroids[i];

                float backY = c.y + PlazaLift + CanopyClearance;
                float frontY = backY - CanopySlope;
                Vector3 bL = new Vector3(c.x - CanopyHalf - CanopyOverhang, backY,  c.z - CanopyHalf - CanopyOverhang);
                Vector3 bR = new Vector3(c.x + CanopyHalf + CanopyOverhang, backY,  c.z - CanopyHalf - CanopyOverhang);
                Vector3 fL = new Vector3(c.x - CanopyHalf - CanopyOverhang, frontY, c.z + CanopyHalf + CanopyOverhang);
                Vector3 fR = new Vector3(c.x + CanopyHalf + CanopyOverhang, frontY, c.z + CanopyHalf + CanopyOverhang);
                MeshPrimitives.AddQuad(verts, tris, bL, bR, fR, fL);  // top
                MeshPrimitives.AddQuad(verts, tris, fL, fR, bR, bL);  // bottom (double-sided)
            }

            return MeshPrimitives.Finalize(verts, tris);
        }
    }
}
