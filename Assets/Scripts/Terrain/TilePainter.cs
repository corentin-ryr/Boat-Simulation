using UnityEngine;

namespace TerrainGrid
{
    // Click-paint primal tiles from gameplay. Casts a ray from the camera through the cursor,
    // hits the terrain MeshCollider, resolves the hit point to (ChunkCoord, polygonIndex) via
    // nearest-centroid lookup on the published render mirror, then enqueues a tile mutation
    // on the streamer. Worker thread applies it next pass; the chunk's Version bumps, the
    // mesh rebuilds, vertex-colour palette flips for that cell — and the paint shows up.
    //
    // Drop this MonoBehaviour onto any GameObject in the scene (the player, the
    // ChunkManager's GameObject, or an empty), assign a ChunkManager reference, and play.
    public class TilePainter : MonoBehaviour
    {
        [Tooltip("Source of the streamer + render mirror. Required.")]
        public ChunkManager chunkManager;

        [Tooltip("Ray origin. If null, falls back to Camera.main at Awake.")]
        public Camera cam;

        [Tooltip("Restrict the raycast to terrain colliders. Leave at -1 (Everything) if your " +
                 "terrain isn't on a dedicated layer.")]
        public LayerMask terrainMask = ~0;

        [Tooltip("Maximum cast distance in world units.")]
        public float maxDistance = 500f;

        [Tooltip("Mouse button (0=left, 1=right, 2=middle) that triggers a paint.")]
        public int paintMouseButton = 0;

        [Tooltip("Key that cycles through TileKind values for the next paint.")]
        public KeyCode cycleKey = KeyCode.T;

        [Tooltip("Show the currently-selected TileKind in the top-left corner.")]
        public bool drawHud = true;

        TileKind current = TileKind.Road;

        void Awake()
        {
            if (cam == null) cam = Camera.main;
        }

        void Update()
        {
            if (chunkManager == null || cam == null) return;

            if (Input.GetKeyDown(cycleKey)) CycleKind();
            if (Input.GetMouseButtonDown(paintMouseButton)) TryPaintUnderCursor();
        }

        void CycleKind()
        {
            // Wrap through every defined TileKind value (Default, Road, Building, ...). Using
            // GetValues keeps this honest if the enum grows.
            var values = System.Enum.GetValues(typeof(TileKind));
            int idx = System.Array.IndexOf(values, current);
            idx = (idx + 1) % values.Length;
            current = (TileKind)values.GetValue(idx);
        }

        void TryPaintUnderCursor()
        {
            Ray ray = cam.ScreenPointToRay(Input.mousePosition);
            if (!Physics.Raycast(ray, out RaycastHit hit, maxDistance, terrainMask)) return;

            // Resolve the hit world position to the chunk it lives in. ChunkManager owns the
            // hex math constants; reading them keeps the painter's bookkeeping aligned with
            // the streamer's coord system without duplicating constants.
            float hexRadius = chunkManager.HexRadius;
            int chunkGridSize = chunkManager.ChunkGridSize;
            ChunkCoord coord = ChunkCoord.FromWorldPos(hit.point, hexRadius, chunkGridSize);

            if (!chunkManager.TryGetPublishedChunk(coord, out PrimalChunk chunk)) return;
            if (chunk.PrimalCellCentroids == null || chunk.PrimalCellCentroids.Length == 0) return;

            // Nearest-centroid lookup. Tiles are small (~hexRadius worth of XZ extent) and
            // chunks contain ~50 cells, so this is a few microseconds — no spatial index
            // needed. Y is ignored because elevation can put the centroid above or below the
            // hit point; the painted tile is the one whose XZ footprint owns the cursor.
            int bestIdx = -1;
            float bestSq = float.PositiveInfinity;
            for (int i = 0; i < chunk.PrimalCellCentroids.Length; i++)
            {
                Vector3 c = chunk.PrimalCellCentroids[i];
                float dx = c.x - hit.point.x;
                float dz = c.z - hit.point.z;
                float sq = dx * dx + dz * dz;
                if (sq < bestSq) { bestSq = sq; bestIdx = i; }
            }
            if (bestIdx < 0) return;

            chunkManager.Streamer.EnqueueTileMutation(
                coord, bestIdx, new TileProperty { Kind = current });
        }

        void OnGUI()
        {
            if (!drawHud) return;
            GUI.color = Color.white;
            GUI.Box(new Rect(10, 10, 220, 26), $"Painting: {current}   ({cycleKey} to cycle)");
        }
    }
}
