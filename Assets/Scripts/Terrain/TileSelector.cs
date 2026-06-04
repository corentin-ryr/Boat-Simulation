using UnityEngine;

namespace TerrainGrid
{
    // RTS-mode tile picker + hover highlight. Reads the cursor via TilePicker each
    // frame and pushes the picked (ChunkCoord, cellIdx) into three global shader
    // uniforms that TerrainElevation.shader reads to tint the hovered cell
    // pixel-exact, with no geometry at all on the C# side.
    //
    // The previous version of this script built a per-frame "hover Mesh" that
    // tried to mimic the terrain's wedge decomposition; it kept sinking below
    // unpainted tiles because the terrain switches between fan triangulation
    // (Default cells) and wedge triangulation (carved cells), and there's no
    // single mesh that exactly matches both. Driving the highlight from the
    // shader instead — using the same nearest-corner barycentric decision the
    // shader already does for tile tinting — eliminates the geometric mismatch
    // entirely: the tint follows the cell's actual rendered pixels.
    //
    // Disabled in gameplay mode; GameModeController toggles the GameObject on
    // Tab. While disabled the script restores _HoverPrimalIdx = -1 so no stale
    // tint survives a mode flip.
    public class TileSelector : MonoBehaviour
    {
        [Tooltip("Source of the streamer + render mirror. Required.")]
        public ChunkManager chunkManager;

        [Tooltip("Camera used for the cursor raycast. Set automatically by " +
                 "GameModeController to the active RTS camera on every mode flip; " +
                 "Inspector value is overwritten.")]
        public Camera cam;

        [Tooltip("Palette that holds the currently-selected TileKind. Required — " +
                 "clicks are dropped without it.")]
        public TilePalette palette;

        [Tooltip("Restrict the raycast to terrain colliders. Leave at -1 (Everything).")]
        public LayerMask terrainMask = ~0;

        [Tooltip("Maximum cast distance in world units.")]
        public float maxDistance = 1000f;

        [Tooltip("Mouse button index that triggers a paint (0 = left).")]
        public int paintMouseButton = 0;

        // Shader-side hover state. -1 means "no cell hovered" and the shader skips
        // the tint entirely. (q, r) is the hovered chunk's axial coord, used by
        // the shader to compare against each chunk's per-renderer _ChunkQR so the
        // same local index doesn't tint cells in every chunk at once.
        static readonly int HoverPrimalIdxId = Shader.PropertyToID("_HoverPrimalIdx");
        static readonly int HoverChunkQRId   = Shader.PropertyToID("_HoverChunkQR");

        // Cached last-picked cell — used to skip redundant SetGlobal calls when
        // the cursor moves within the same cell (the dominant case).
        ChunkCoord lastCoord;
        int lastCellIdx = -1;
        bool hasLastCell;

        void OnEnable()
        {
            // Reset cache on entry so the next pick definitely pushes.
            hasLastCell = false;
            lastCellIdx = -1;
        }

        void OnDisable()
        {
            // Clear the hover so the tint doesn't survive into gameplay mode.
            Shader.SetGlobalFloat(HoverPrimalIdxId, -1f);
            hasLastCell = false;
            lastCellIdx = -1;
        }

        void Update()
        {
            if (chunkManager == null || cam == null) return;

            bool hit = TilePicker.PickUnderCursor(
                cam, chunkManager, terrainMask, maxDistance,
                out ChunkCoord coord, out int cellIdx, out Vector3 _);

            if (!hit)
            {
                if (hasLastCell)
                {
                    Shader.SetGlobalFloat(HoverPrimalIdxId, -1f);
                    hasLastCell = false;
                    lastCellIdx = -1;
                }
                return;
            }

            bool sameCell = hasLastCell && coord == lastCoord && cellIdx == lastCellIdx;
            if (!sameCell)
            {
                Shader.SetGlobalFloat(HoverPrimalIdxId, cellIdx);
                Shader.SetGlobalVector(HoverChunkQRId, new Vector4(coord.q, coord.r, 0f, 0f));
                lastCoord = coord;
                lastCellIdx = cellIdx;
                hasLastCell = true;
            }

            // Paint: left-click enqueues a mutation with the palette's current kind.
            // Skip if the cursor is over the palette HUD itself — otherwise clicking
            // a palette button paints the underlying tile too.
            //
            // TilePlacement.IsValid filters out illegal targets BEFORE the mutation
            // reaches the streamer — Lighthouse/Dock require a coastal cell or
            // adjacency to an existing Lighthouse/Dock. Invalid clicks are silently
            // dropped; v1 has no on-screen feedback yet, but the outline already
            // shows the cell being targeted so the user can move to a valid spot.
            if (Input.GetMouseButtonDown(paintMouseButton)
                && palette != null
                && !palette.IsMouseOverPalette()
                && TilePlacement.IsValid(palette.Current, coord, cellIdx, chunkManager))
            {
                chunkManager.Streamer.EnqueueTileMutation(
                    coord, cellIdx, new TileProperty { Kind = palette.Current });
            }
        }
    }
}
