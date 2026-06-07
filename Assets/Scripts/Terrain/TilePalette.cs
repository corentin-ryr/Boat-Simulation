using UnityEngine;

namespace TerrainGrid
{
    // Owns the currently-selected TileKind for the painter system. Shared by
    // TilePainter (gameplay-mode left-click paint) and TileSelector (RTS-mode
    // click-to-paint). Decoupling the "what am I painting?" state from either
    // picker means flipping between modes preserves the selection, and there's a
    // single source of truth for the HUD to render.
    //
    // The HUD + hotkeys are driven by a hardcoded `Entries` table at the top
    // of the class. Adding a new TileKind = add one entry here + add its
    // overlay in ChunkSurface. No ScriptableObjects, no Inspector wiring
    // beyond the per-kind Materials on ChunkManager.
    //
    // Controls:
    //   T              cycle through all entries in table order
    //   1..7           jump directly to the kind whose hotkey matches
    //   GUI button     same effect as the hotkey
    public class TilePalette : MonoBehaviour
    {
        [Tooltip("Initial selection at play start.")]
        public TileKind initialKind = TileKind.Road;

        [Tooltip("Show the palette bar in the top-left corner. Off → palette still " +
                 "drives the selection via keyboard, just no on-screen UI.")]
        public bool drawHud = true;

        [Tooltip("Key that cycles through the table in order.")]
        public KeyCode cycleKey = KeyCode.T;

        // Read by TilePainter and TileSelector when emitting an EnqueueTileMutation.
        public TileKind Current { get; private set; } = TileKind.Default;

        // Hardcoded palette table — adding a TileKind = one new entry here.
        // Order is the on-screen button order and the cycle order. Display
        // name is shown on the button; hotkey is the keyboard shortcut.
        struct Entry
        {
            public TileKind kind;
            public string name;
            public KeyCode hotkey;
            public Entry(TileKind k, string n, KeyCode h) { kind = k; name = n; hotkey = h; }
        }
        static readonly Entry[] Entries =
        {
            new Entry(TileKind.Default,    "Default",    KeyCode.Alpha1),
            new Entry(TileKind.Road,       "Road",       KeyCode.Alpha2),
            new Entry(TileKind.Building,   "Building",   KeyCode.Alpha3),
            new Entry(TileKind.Field,      "Field",      KeyCode.Alpha4),
            new Entry(TileKind.Market,     "Market",     KeyCode.Alpha5),
            new Entry(TileKind.Lighthouse, "Lighthouse", KeyCode.Alpha6),
            new Entry(TileKind.Dock,       "Dock",       KeyCode.Alpha7),
        };

        // Layout constants shared between OnGUI drawing and the screen-space hit-test.
        const float BtnW   = 100f;
        const float BtnH   = 26f;
        const float Pad    = 4f;
        const float Margin = 10f;
        static readonly float PaletteWidth  = Entries.Length * BtnW + (Entries.Length - 1) * Pad;
        static readonly float PaletteHeight = BtnH + 22f;

        // True while the cursor sits over the palette HUD rect. Used by TileSelector
        // to suppress paint clicks that the OnGUI buttons should own. Computed in
        // Input-coords (Y origin = bottom of screen) regardless of whether the
        // GUI is drawn this frame.
        public bool IsMouseOverPalette()
        {
            if (!drawHud) return false;
            float mx = Input.mousePosition.x;
            float my = Input.mousePosition.y;
            float yTop = Screen.height - Margin;
            float yBot = yTop - PaletteHeight;
            return mx >= Margin && mx <= Margin + PaletteWidth
                && my >= yBot   && my <= yTop;
        }

        void Awake()
        {
            Current = initialKind;
        }

        void Update()
        {
            if (Input.GetKeyDown(cycleKey)) CycleKind();
            for (int i = 0; i < Entries.Length; i++)
            {
                if (Input.GetKeyDown(Entries[i].hotkey)) { Current = Entries[i].kind; break; }
            }
        }

        void CycleKind()
        {
            int curIdx = 0;
            for (int i = 0; i < Entries.Length; i++)
                if (Entries[i].kind == Current) { curIdx = i; break; }
            int next = (curIdx + 1) % Entries.Length;
            Current = Entries[next].kind;
        }

        void OnGUI()
        {
            if (!drawHud) return;

            float x = Margin;
            float y = Margin;
            for (int i = 0; i < Entries.Length; i++)
            {
                DrawKindButton(new Rect(x, y, BtnW, BtnH), $"{i + 1}: {Entries[i].name}", Entries[i].kind);
                x += BtnW + Pad;
            }
            GUI.color = new Color(1, 1, 1, 0.7f);
            GUI.Label(new Rect(Margin, y + BtnH + 2f, 400f, 18f),
                      $"Press {cycleKey} to cycle / 1-{Entries.Length} to jump");
            GUI.color = Color.white;
        }

        void DrawKindButton(Rect rect, string label, TileKind kind)
        {
            bool selected = kind == Current;
            Color prev = GUI.color;
            GUI.color = selected ? new Color(1f, 1f, 0.4f, 1f) : Color.white;
            if (GUI.Button(rect, label)) Current = kind;
            GUI.color = prev;
        }
    }
}
