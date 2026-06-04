using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.Presenters
{
    // Dock: a flat plank-board raised slightly above the terrain, with four
    // short support posts hanging below at the corners. The board is wider on
    // X than Z, giving an implicit "long axis pointing into the water"
    // orientation once we add per-cell coastal-facing rotation.
    [CreateAssetMenu(menuName = "Terrain/Presenter/Dock", fileName = "DockPresenter")]
    public class DockPresenter : TilePresenter
    {
        public float halfX = 0.45f;
        public float halfZ = 0.30f;
        public float boardThickness = 0.08f;
        public float lift = 0.06f;
        public float postHalf = 0.05f;
        public float postDrop = 0.30f;

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

                // Board.
                Vector3 boardMin = c + new Vector3(-halfX, lift, -halfZ);
                Vector3 boardMax = c + new Vector3(halfX, lift + boardThickness, halfZ);
                MeshUtils.AddBox(v, t, boardMin, boardMax);

                // Support posts at each corner.
                float postTop = c.y + lift;
                float postBot = c.y + lift - postDrop;
                Vector2[] corners = {
                    new Vector2(-halfX + postHalf, -halfZ + postHalf),
                    new Vector2(halfX - postHalf, -halfZ + postHalf),
                    new Vector2(-halfX + postHalf, halfZ - postHalf),
                    new Vector2(halfX - postHalf, halfZ - postHalf),
                };
                foreach (Vector2 cor in corners)
                {
                    Vector3 pMin = new Vector3(c.x + cor.x - postHalf, postBot, c.z + cor.y - postHalf);
                    Vector3 pMax = new Vector3(c.x + cor.x + postHalf, postTop, c.z + cor.y + postHalf);
                    MeshUtils.AddBox(v, t, pMin, pMax);
                }
            }

            return MeshUtils.Finalize(v, t);
        }
    }
}
