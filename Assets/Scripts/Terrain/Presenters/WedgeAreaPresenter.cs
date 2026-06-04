using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Presenters
{
    // Flat-area presenter: covers each cell of this kind with a wedge-decomposed
    // quad sheet at terrain Y (+ optional lift to escape z-fighting). This is the
    // shape Field uses today, and it's the right reuse target for any future
    // flat-ground kinds (Plaza, Sand, Pavement) — clone the asset, tweak the lift
    // and material, drop into the catalog. Zero new C#.
    //
    // The terrain mesh side carves the wedge of every cell whose kind has
    // `TileDefinition.carvesTerrain = true`, so this overlay's quad IS the
    // visible cell — not a decal on top. The tiny `lift` value is just to win
    // the z-fight at the carved boundary; bigger lifts will make the area
    // visibly float above the terrain, which is fine for stylised gameplay
    // visuals if you want it.
    [CreateAssetMenu(menuName = "Terrain/Presenter/Wedge Area", fileName = "WedgeAreaPresenter")]
    public class WedgeAreaPresenter : TilePresenter
    {
        [Tooltip("Metres the wedge is lifted above the terrain centroid to avoid z-fighting " +
                 "with the carved hole's edge. 0.03 is fine at hexRadius=1.")]
        public float lift = 0.03f;

        protected override Mesh BuildMesh(in PresenterContext ctx)
        {
            var vertices = new List<Vector3>();
            var triangles = new List<int>();

            PrimalChunk primal = ctx.primal;
            List<Polygon> polygons = primal.Dual;
            int[] dvpi = primal.DualVertexPrimalIndex;
            TileProperty[] tp = primal.TileProperties;
            if (polygons == null || dvpi == null || tp == null) return MeshUtils.EmptyMesh();

            Vector3 liftV = new Vector3(0f, lift, 0f);
            TileKind kind = ctx.kind;

            int cornerCursor = 0;
            foreach (Polygon poly in polygons)
            {
                Vector3[] verts = poly.GetVerticesPosition();
                int n = verts.Length;
                if (n < 3) { cornerCursor += n; continue; }

                bool[] isKind = new bool[n];
                bool any = false;
                for (int i = 0; i < n; i++)
                {
                    int pi = cornerCursor + i < dvpi.Length ? dvpi[cornerCursor + i] : -1;
                    bool k = pi >= 0 && pi < tp.Length && tp[pi].Kind == kind;
                    isKind[i] = k;
                    if (k) any = true;
                }
                cornerCursor += n;
                if (!any) continue;

                Vector3 V = Vector3.zero;
                for (int i = 0; i < n; i++) V += verts[i];
                V /= n;

                for (int i = 0; i < n; i++)
                {
                    if (!isKind[i]) continue;
                    int iPrev = (i + n - 1) % n;
                    int iNext = (i + 1) % n;
                    Vector3 Mprev = 0.5f * (verts[iPrev] + verts[i]);
                    Vector3 Mnext = 0.5f * (verts[i] + verts[iNext]);

                    // Same (c_i, V, M_next) / (c_i, M_prev, V) winding as the
                    // terrain mesh — guarantees the overlay normal points up
                    // and matches the terrain orientation.
                    int b = vertices.Count;
                    vertices.Add(verts[i] + liftV);
                    vertices.Add(Mnext + liftV);
                    vertices.Add(V + liftV);
                    vertices.Add(Mprev + liftV);
                    triangles.Add(b + 0); triangles.Add(b + 2); triangles.Add(b + 1);
                    triangles.Add(b + 0); triangles.Add(b + 3); triangles.Add(b + 2);
                }
            }

            return MeshUtils.Finalize(vertices, triangles);
        }
    }
}
