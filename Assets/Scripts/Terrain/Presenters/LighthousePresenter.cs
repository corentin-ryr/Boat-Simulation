using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Presenters
{
    // Lighthouse: three stacked tapered prisms (wide → mid → narrow) plus a
    // small cap cone, around the cell's primal centroid. Tall enough to read
    // from any RTS-camera zoom level. Standalone — the terrain wedge under
    // the base stays visible.
    [CreateAssetMenu(menuName = "Terrain/Presenter/Lighthouse", fileName = "LighthousePresenter")]
    public class LighthousePresenter : TilePresenter
    {
        public float baseRadius = 0.30f;
        public float midRadius = 0.22f;
        public float topRadius = 0.16f;
        public float baseHeight = 1.6f;
        public float midHeight = 1.6f;
        public float topHeight = 1.2f;
        public float capHeight = 0.35f;
        public int sides = 8;

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

                float y0 = c.y;
                float y1 = y0 + baseHeight;
                float y2 = y1 + midHeight;
                float y3 = y2 + topHeight;
                float y4 = y3 + capHeight;

                MeshUtils.AddTaperedPrism(v, t, c.x, c.z, y0, baseRadius, y1, baseRadius, sides);
                MeshUtils.AddTaperedPrism(v, t, c.x, c.z, y1, baseRadius, y2, midRadius, sides);
                MeshUtils.AddTaperedPrism(v, t, c.x, c.z, y2, midRadius, y3, topRadius, sides);
                MeshUtils.AddTaperedPrism(v, t, c.x, c.z, y3, topRadius, y4, 0.02f, sides);
            }

            return MeshUtils.Finalize(v, t);
        }
    }
}
