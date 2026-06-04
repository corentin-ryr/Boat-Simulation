using UnityEngine;

namespace TerrainGrid
{
    // Click-paint primal tiles from the GAMEPLAY camera. Casts a ray through the
    // cursor via the shared TilePicker, then enqueues a tile mutation using the
    // currently-selected TileKind from the shared TilePalette.
    //
    // The actual picker logic (raycast → ChunkCoord → nearest-centroid) lives in
    // TilePicker, so the gameplay and RTS painters walk the exact same code path
    // — no behavioural drift between modes. The selected TileKind lives on
    // TilePalette, so flipping between modes preserves the selection and there's
    // a single source of truth for what gets painted.
    //
    // GameModeController enables this MonoBehaviour while in gameplay mode and
    // disables it in RTS mode (where TileSelector owns the painter role).
    public class TilePainter : MonoBehaviour
    {
        [Tooltip("Source of the streamer + render mirror. Required.")]
        public ChunkManager chunkManager;

        [Tooltip("Ray origin. If null, falls back to Camera.main at Awake. Should be " +
                 "the GAMEPLAY camera while this script is enabled.")]
        public Camera cam;

        [Tooltip("Source of truth for the currently-selected TileKind. Required.")]
        public TilePalette palette;

        [Tooltip("Restrict the raycast to terrain colliders. Leave at -1 (Everything) if " +
                 "your terrain isn't on a dedicated layer.")]
        public LayerMask terrainMask = ~0;

        [Tooltip("Maximum cast distance in world units.")]
        public float maxDistance = 500f;

        [Tooltip("Mouse button (0=left, 1=right, 2=middle) that triggers a paint.")]
        public int paintMouseButton = 0;

        void Awake()
        {
            if (cam == null) cam = Camera.main;
        }

        void Update()
        {
            if (chunkManager == null || cam == null || palette == null) return;
            if (!Input.GetMouseButtonDown(paintMouseButton)) return;

            if (!TilePicker.PickUnderCursor(
                    cam, chunkManager, terrainMask, maxDistance,
                    out ChunkCoord coord, out int cellIdx, out Vector3 _))
                return;

            // Drop illegal placements before enqueueing (Lighthouse/Dock require coast).
            // Gameplay mode has no on-screen feedback for the rejection; the click just
            // does nothing, same as missing the terrain.
            if (!TilePlacement.IsValid(palette.Current, coord, cellIdx, chunkManager))
                return;

            chunkManager.Streamer.EnqueueTileMutation(
                coord, cellIdx, new TileProperty { Kind = palette.Current });
        }
    }
}
