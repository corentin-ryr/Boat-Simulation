using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Presenters
{
    // Raised-roof building overlay: each Building cell's wedge is lifted to a
    // per-cell roof height (`primalCentroid.y + height`), and vertical walls
    // descend from the roof to the terrain (or to the lower-roofed neighbour's
    // roof, for adjacent buildings with stepped heights). Adjacent Building
    // cells with the same roofY merge into one flat roof with no internal wall.
    //
    // RequiresCollider=true — terrain is carved away under the footprint and the
    // walls need to block the player from walking through.
    [CreateAssetMenu(menuName = "Terrain/Presenter/Building", fileName = "BuildingPresenter")]
    public class BuildingPresenter : TilePresenter
    {
        [Tooltip("Roof height above the cell's primal centroid Y. Per-cell flat — adjacent " +
                 "Building cells on a slope terrace vertically at the cell boundary, with " +
                 "the higher cell emitting the step wall.")]
        public float height = 3.0f;

        public override bool RequiresCollider => true;

        protected override Mesh BuildMesh(in PresenterContext ctx)
        {
            var verts = new List<Vector3>();
            var tris = new List<int>();

            PrimalChunk primal = ctx.primal;
            List<Polygon> dualPolys = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            Vector3[] centroids = primal.PrimalCellCentroids;
            if (dualPolys == null || tp == null || centroids == null) return MeshUtils.EmptyMesh();

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
                    if (b) { roofY[i] = centroids[pi].y + height; anyBuilding = true; }
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

                    MeshUtils.AddTriangle(verts, tris, ciHi, Vhi, MnextHi);
                    MeshUtils.AddTriangle(verts, tris, ciHi, MprevHi, Vhi);

                    // V→Mprev wall. Neighbour wedge iPrev controls the bottom Y
                    // along the segment: Building neighbour → its flat roofY;
                    // non-Building → per-endpoint terrain Y so the wall tracks the
                    // slope at its base.
                    {
                        bool nbBuilding = isBuilding[iPrev];
                        float bottomAtV = nbBuilding ? roofY[iPrev] : V.y;
                        float bottomAtM = nbBuilding ? roofY[iPrev] : Mprev.y;
                        if (ry > bottomAtV || ry > bottomAtM)
                        {
                            Vector3 VlowP = new Vector3(V.x, bottomAtV, V.z);
                            Vector3 MprevLow = new Vector3(Mprev.x, bottomAtM, Mprev.z);
                            EmitWall(verts, tris, VlowP, MprevLow, MprevHi, Vhi, ci);
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
                            EmitWall(verts, tris, MnextLow, VlowN, Vhi, MnextHi, ci);
                        }
                    }
                }
            }

            return MeshUtils.Finalize(verts, tris);
        }

        // Vertical wall facing outward (away from `inwardRef`, the wedge's
        // road-side corner). Single-sided — the inside is carved-terrain void
        // with nothing to render against, so the back face is invisible by
        // construction. Picks winding from the cross-product sign against the
        // outward direction.
        static void EmitWall(List<Vector3> verts, List<int> tris,
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
                MeshUtils.AddQuad(verts, tris, aLow, bLow, bHigh, aHigh);
            else
                MeshUtils.AddQuad(verts, tris, bLow, aLow, aHigh, bHigh);
        }
    }
}
