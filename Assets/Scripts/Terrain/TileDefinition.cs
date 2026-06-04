using UnityEngine;

namespace TerrainGrid
{
    // Static, authored metadata for one TileKind. Everything the rest of the
    // codebase needs to know about a kind — display name, hotkey, palette
    // colour, material, how to render it, where it's allowed — lives here as
    // a single .asset file. Adding a new TileKind means: extend the enum byte
    // (for snapshot stability), make a new TileDefinition asset, drop it into
    // the TileCatalog. No C# code in the common case.
    //
    // Two slots can be null:
    //   - `presenter == null` → kind has no visual overlay (Default). The
    //     terrain mesh paints the cell directly.
    //   - `placementRule == null` → no restrictions; the user can paint this
    //     kind anywhere.
    [CreateAssetMenu(menuName = "Terrain/Tile Definition", fileName = "Tile")]
    public class TileDefinition : ScriptableObject
    {
        [Tooltip("Enum byte this definition describes. Must be unique across the catalog. " +
                 "Stable identifier — survives save/load via TileProperty.Kind.")]
        public TileKind kind;

        [Tooltip("Shown on the palette button and (eventually) info HUDs.")]
        public string displayName;

        [Tooltip("Keyboard key that jumps the palette directly to this kind. Use 1-9 or " +
                 "letters; the palette's cycle key (default T) walks through every kind in " +
                 "catalog order.")]
        public KeyCode hotkey = KeyCode.None;

        [Tooltip("Highlight tint applied to the palette button when this kind is selected. " +
                 "Used in IMGUI; no shader-side meaning yet.")]
        public Color paletteSwatch = new Color(1f, 1f, 0.4f, 1f);

        [Tooltip("Material applied to the presenter's overlay MeshRenderer. Null is legal " +
                 "for kinds with no presenter (Default).")]
        public Material material;

        [Tooltip("How the cell renders. Null = no overlay; the terrain mesh handles it. " +
                 "Pluggable: the same WedgeAreaPresenter asset can be reused for any new " +
                 "flat-area tile (Sand, Plaza) by cloning + tweaking parameters.")]
        public TilePresenter presenter;

        [Tooltip("Where the user can paint this kind. Null = anywhere. Same rule asset can " +
                 "be reused across kinds (e.g. one CoastalOrAdjacentRule covers both " +
                 "Lighthouse and Dock).")]
        public TilePlacementRule placementRule;

        [Tooltip("True if the terrain mesh should SKIP the wedge of cells of this kind — " +
                 "the presenter's mesh paints the cell pixel-for-pixel. Use for area " +
                 "presenters (Road, Field, Building) that cover the natural ground. False " +
                 "for point-structure presenters (Market, Lighthouse, Dock) that sit on " +
                 "top of the visible terrain.")]
        public bool carvesTerrain;
    }
}
