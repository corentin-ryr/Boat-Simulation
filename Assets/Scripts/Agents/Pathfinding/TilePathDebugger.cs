using TerrainGrid;
using UnityEngine;

namespace Agents.Pathfinding
{
    // Scene-time visualizer for TilePathfinder. Click a tile in any input
    // mode to set the "from" anchor; Shift-click to set "to" and run the
    // pathfinder. The resulting path draws every frame via Debug.DrawLine
    // (visible in the Scene view; in the Game view enable Gizmos).
    //
    // Wired through TilePicker so the cell-resolution math stays in
    // exactly one place. No physical raycast layer required beyond what
    // TilePainter / TileSelector already use.
    public class TilePathDebugger : MonoBehaviour
    {
        [Tooltip("ChunkManager whose mirror this debugger reads from.")]
        public ChunkManager chunkManager;

        [Tooltip("Camera used for the cursor raycast. Defaults to Camera.main.")]
        public Camera cam;

        [Tooltip("Layer mask for the terrain hit. Match the painters' setup.")]
        public LayerMask terrainMask = ~0;

        [Tooltip("Maximum raycast distance.")]
        public float maxDistance = 500f;

        [Tooltip("Draw the path in this colour. Bright green by default so it pops " +
                 "over the terrain tint.")]
        public Color pathColour = new Color(0.2f, 1f, 0.3f);

        [Tooltip("Lift the drawn line by this much in world Y so it doesn't z-fight " +
                 "with the wedge tops.")]
        public float drawLift = 0.5f;

        [Tooltip("Off by default so the debugger doesn't grab clicks meant for the " +
                 "painters. F6 hot-key toggles at runtime; when on, the next left-click " +
                 "sets 'from' and shift-left-click runs the pathfinder.")]
        public bool enableInput = false;
        public KeyCode toggleKey = KeyCode.F6;

        bool hasFrom;
        ChunkCoord fromCoord;
        int fromCell;

        TilePath lastPath;
        string lastStatus = "—";

        TilePathCache cache;

        void Awake()
        {
            if (chunkManager == null) chunkManager = FindObjectOfType<ChunkManager>();
            if (cam == null) cam = Camera.main;
        }

        void Start()
        {
            if (chunkManager != null) cache = new TilePathCache(chunkManager, capacity: 64);
        }

        void Update()
        {
            if (Input.GetKeyDown(toggleKey)) enableInput = !enableInput;
            if (!enableInput) return;
            if (cam == null || chunkManager == null) return;

            // Left-click without Shift = from anchor.
            // Left-click with Shift     = to anchor (runs the pathfinder).
            if (Input.GetMouseButtonDown(0))
            {
                bool shift = Input.GetKey(KeyCode.LeftShift) || Input.GetKey(KeyCode.RightShift);
                if (!TilePicker.PickUnderCursor(cam, chunkManager, terrainMask, maxDistance,
                                                out ChunkCoord c, out int cell, out _))
                    return;
                if (!shift)
                {
                    fromCoord = c; fromCell = cell; hasFrom = true;
                    lastPath = null;
                    lastStatus = $"from = {c}:{cell}";
                }
                else if (hasFrom)
                {
                    if (cache == null) cache = new TilePathCache(chunkManager, capacity: 64);
                    lastPath = cache.GetOrFind(fromCoord, fromCell, c, cell);
                    lastStatus = lastPath != null
                        ? $"path: {lastPath.Length} steps, cost {lastPath.TotalCost:0.0}"
                        : "path: NONE (unreachable / out of nodes)";
                }
            }

            DrawPath();
        }

        void DrawPath()
        {
            if (lastPath == null || lastPath.Length < 2) return;
            Vector3 lift = new Vector3(0f, drawLift, 0f);
            for (int i = 0; i < lastPath.Length - 1; i++)
            {
                Vector3 a = lastPath[i].WorldPos + lift;
                Vector3 b = lastPath[i + 1].WorldPos + lift;
                Debug.DrawLine(a, b, pathColour);
            }
        }

        void OnGUI()
        {
            const float w = 280f;
            float x = 10f;
            float y = Screen.height - 50f;
            GUI.Box(new Rect(x, y, w, 40f), "");
            GUI.Label(new Rect(x + 6f, y + 2f, w - 12f, 18f),
                      $"PathDebugger  [{toggleKey} toggle]  enabled={enableInput}");
            GUI.Label(new Rect(x + 6f, y + 20f, w - 12f, 18f), lastStatus);
        }
    }
}
