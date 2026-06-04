using UnityEngine;

namespace TerrainGrid
{
    // Owns the currently-selected TileKind for the painter system. Shared by
    // TilePainter (gameplay-mode left-click paint) and TileSelector (RTS-mode
    // click-to-paint). Decoupling the "what am I painting?" state from either
    // picker means flipping between modes preserves the selection, and there's a
    // single source of truth for the HUD to render.
    //
    // The HUD + hotkeys are driven entirely by the TileCatalog asset on
    // ChunkManager — no hardcoded TileKind references here. Adding a new tile
    // appends one button + one hotkey to the bar without code edits.
    //
    // Controls:
    //   T              cycle through every kind in catalog order
    //   def.hotkey     jump directly to the kind whose definition has it set
    //   GUI button     same effect as the hotkey
    //
    // The HUD is OnGUI for parity with the rest of the IMGUI surface in this
    // project (on-screen profiler stats overlay, gameplay TilePainter HUD).
    public class TilePalette : MonoBehaviour
    {
        [Tooltip("ChunkManager to read the TileCatalog from. Set automatically by " +
                 "PlayerSpawner / GameModeController if left null. Without a catalog, the " +
                 "HUD is empty and the palette serves Current = TileKind.Default.")]
        public ChunkManager chunkManager;

        [Tooltip("Initial selection at play start. Must match an entry in the catalog.")]
        public TileKind initialKind = TileKind.Road;

        [Tooltip("Show the palette bar in the top-left corner. Off → palette still " +
                 "drives the selection via keyboard, just no on-screen UI.")]
        public bool drawHud = true;

        [Tooltip("Key that cycles through all catalog entries (legacy behaviour).")]
        public KeyCode cycleKey = KeyCode.T;

        // Read by TilePainter and TileSelector when emitting an EnqueueTileMutation.
        public TileKind Current { get; private set; } = TileKind.Default;

        // Layout constants shared between OnGUI drawing and the screen-space hit-test.
        const float BtnW = 100f;
        const float BtnH = 26f;
        const float Pad  = 4f;
        const float Margin = 10f;

        TileCatalog Catalog => chunkManager != null ? chunkManager.Catalog : null;

        // True while the cursor sits over the palette HUD rect. Used by TileSelector
        // to suppress paint clicks that the OnGUI buttons should own. Computed in
        // Input-coords (Y origin = bottom of screen) regardless of whether the
        // GUI is drawn this frame.
        public bool IsMouseOverPalette()
        {
            if (!drawHud) return false;
            TileCatalog cat = Catalog;
            int count = cat != null ? cat.Definitions.Count : 0;
            if (count == 0) return false;

            float w = count * BtnW + (count - 1) * Pad;
            float h = BtnH + 22f;
            float mx = Input.mousePosition.x;
            float my = Input.mousePosition.y;
            float yTop = Screen.height - Margin;
            float yBot = yTop - h;
            return mx >= Margin && mx <= Margin + w
                && my >= yBot && my <= yTop;
        }

        void Awake()
        {
            Current = initialKind;
            // Scene compatibility: when the field is left blank on an existing scene,
            // grab the first ChunkManager in the scene. Single-ChunkManager scenes (the
            // current setup) work without explicit Inspector wiring.
            if (chunkManager == null)
#if UNITY_2023_1_OR_NEWER
                chunkManager = Object.FindFirstObjectByType<ChunkManager>();
#else
                chunkManager = Object.FindObjectOfType<ChunkManager>();
#endif
        }

        void Update()
        {
            TileCatalog cat = Catalog;
            if (cat == null) return;

            if (Input.GetKeyDown(cycleKey)) CycleKind(cat);

            // Per-definition hotkey dispatch. KeyCode.None entries are silently
            // ignored; multiple definitions with the same hotkey resolve to the
            // first match in catalog order.
            foreach (TileDefinition def in cat.Definitions)
            {
                if (def == null || def.hotkey == KeyCode.None) continue;
                if (Input.GetKeyDown(def.hotkey)) { Current = def.kind; break; }
            }
        }

        void CycleKind(TileCatalog cat)
        {
            int n = cat.Definitions.Count;
            if (n == 0) return;
            int curIdx = 0;
            for (int i = 0; i < n; i++)
                if (cat.Definitions[i] != null && cat.Definitions[i].kind == Current) { curIdx = i; break; }
            int next = (curIdx + 1) % n;
            TileDefinition def = cat.Definitions[next];
            if (def != null) Current = def.kind;
        }

        void OnGUI()
        {
            if (!drawHud) return;
            TileCatalog cat = Catalog;
            if (cat == null) return;

            float x = Margin;
            float y = Margin;

            foreach (TileDefinition def in cat.Definitions)
            {
                if (def == null) continue;
                string label = HotkeyLabel(def);
                DrawKindButton(new Rect(x, y, BtnW, BtnH), label, def);
                x += BtnW + Pad;
            }

            GUI.color = new Color(1, 1, 1, 0.7f);
            GUI.Label(new Rect(Margin, y + BtnH + 2f, 400f, 18f),
                      $"Press {cycleKey} to cycle / hotkeys to jump");
            GUI.color = Color.white;
        }

        static string HotkeyLabel(TileDefinition def)
        {
            // KeyCode.None shows as plain name; numeric / alpha keys get a "1: " prefix.
            if (def.hotkey == KeyCode.None) return def.displayName;
            return $"{HotkeyText(def.hotkey)}: {def.displayName}";
        }

        static string HotkeyText(KeyCode k)
        {
            // Strip the "Alpha" prefix Unity uses for the top-row digits so "1" shows
            // instead of "Alpha1". Letter keys come through cleanly already.
            string s = k.ToString();
            if (s.StartsWith("Alpha")) return s.Substring(5);
            return s;
        }

        void DrawKindButton(Rect rect, string label, TileDefinition def)
        {
            bool selected = def.kind == Current;
            Color prev = GUI.color;
            GUI.color = selected ? def.paletteSwatch : Color.white;
            if (GUI.Button(rect, label)) Current = def.kind;
            GUI.color = prev;
        }
    }
}
