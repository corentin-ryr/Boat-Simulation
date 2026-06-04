using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Presenters
{
    // Townscaper-style road overlay: each Road cell's wedge is sunk by `sinkDepth`
    // (so the surface sits below the terrain that's been carved away for it), and
    // any wedge boundary that touches a non-Road wedge of the same dual polygon
    // gets a vertical curb running from the sunk Y up to the terrain Y. Adjacent
    // Road wedges share their midpoint at the same sunken height, so the road
    // surface is seamless across cell joins.
    //
    // RequiresCollider=true because the terrain mesh has a hole where Road cells
    // sit — without this overlay's collider, the player would fall through.
    [CreateAssetMenu(menuName = "Terrain/Presenter/Road", fileName = "RoadPresenter")]
    public class RoadPresenter : TilePresenter
    {
        [Tooltip("Metres the road surface sits below the terrain it carves out. 0.06 reads " +
                 "as a believable street with curbs at hexRadius=1.")]
        public float sinkDepth = 0.06f;

        public override bool RequiresCollider => true;

        protected override Mesh BuildMesh(in PresenterContext ctx)
        {
            var verts = new List<Vector3>();
            var tris = new List<int>();

            PrimalChunk primal = ctx.primal;
            List<Polygon> dualPolys = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            if (dualPolys == null || tp == null) return MeshUtils.EmptyMesh();

            Vector3 sink = new Vector3(0, sinkDepth, 0);

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

                    // Sunk road surface — same winding as the terrain wedge it
                    // replaces, so normals point up.
                    MeshUtils.AddTriangle(verts, tris, ciLow, Vlow, MnextLow);
                    MeshUtils.AddTriangle(verts, tris, ciLow, MprevLow, Vlow);

                    // Curb walls on interior wedge boundaries that face a non-Road
                    // wedge. Double-sided so the curb is visible from either side.
                    if (!isRoad[iPrev]) EmitCurb(verts, tris, Vlow, MprevLow, Mprev, V);
                    if (!isRoad[iNext]) EmitCurb(verts, tris, MnextLow, Vlow, V, Mnext);
                }
            }

            return MeshUtils.Finalize(verts, tris);
        }

        // Vertical curb between sunk (aLow, bLow) and surface (bHigh, aHigh).
        // Double-sided — RTS view and on-foot view both want to see it.
        static void EmitCurb(List<Vector3> verts, List<int> tris,
                             Vector3 aLow, Vector3 bLow, Vector3 bHigh, Vector3 aHigh)
        {
            MeshUtils.AddQuad(verts, tris, aLow, bLow, bHigh, aHigh);
            MeshUtils.AddQuad(verts, tris, bLow, aLow, aHigh, bHigh);
        }
    }
}
