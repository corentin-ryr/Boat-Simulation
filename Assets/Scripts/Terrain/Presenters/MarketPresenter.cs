using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Presenters
{
    // Market stall: a low square base + a slanted thin awning above and
    // overhanging on the front (+Z). Placed at the cell's primal centroid; the
    // terrain wedge underneath stays visible because TileDefinition.carvesTerrain
    // is false for this kind.
    //
    // Each cell of this kind in the chunk produces an independent stall in the
    // combined chunk mesh — no instancing yet, since the mesh is small.
    [CreateAssetMenu(menuName = "Terrain/Presenter/Market", fileName = "MarketPresenter")]
    public class MarketPresenter : TilePresenter
    {
        public float baseHalf = 0.28f;
        public float baseHeight = 0.35f;
        public float awningHeight = 0.55f;
        public float awningOverhang = 0.10f;

        protected override Mesh BuildMesh(in PresenterContext ctx)
        {
            var v = new List<Vector3>();
            var t = new List<int>();
            TileProperty[] tp = ctx.primal.TileProperties;
            Vector3[] centroids = ctx.primal.PrimalCellCentroids;
            if (tp == null || centroids == null) return MeshUtils.EmptyMesh();

            for (int i = 0; i < tp.Length && i < centroids.Length; i++)
            {
                if (tp[i].Kind != ctx.kind) continue;
                Vector3 c = centroids[i];

                // Stall base: axis-aligned box from terrain Y up to baseHeight.
                MeshUtils.AddBox(v, t,
                                 c + new Vector3(-baseHalf, 0f, -baseHalf),
                                 c + new Vector3(baseHalf, baseHeight, baseHalf));

                // Awning: tilted slab above the stall, overhanging on +Z. The
                // height drop from back to front gives the "leaning forward"
                // silhouette that makes the shape read as a market.
                float aBackY = baseHeight + awningHeight;
                float aFrontY = baseHeight + awningHeight * 0.55f;
                Vector3 bL = c + new Vector3(-baseHalf - awningOverhang, aBackY, -baseHalf - awningOverhang);
                Vector3 bR = c + new Vector3(baseHalf + awningOverhang, aBackY, -baseHalf - awningOverhang);
                Vector3 fL = c + new Vector3(-baseHalf - awningOverhang, aFrontY, baseHalf + awningOverhang);
                Vector3 fR = c + new Vector3(baseHalf + awningOverhang, aFrontY, baseHalf + awningOverhang);
                MeshUtils.AddQuad(v, t, bL, bR, fR, fL);    // top face
                MeshUtils.AddQuad(v, t, fL, fR, bR, bL);    // bottom face (so the awning isn't one-sided)
            }

            return MeshUtils.Finalize(v, t);
        }
    }
}
