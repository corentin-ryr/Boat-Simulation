using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Per-kind vertex-colour multiplier table for the terrain shader.
    //
    // The terrain shader (Assets/Shader/TerrainElevation.shader) uses
    // TEXCOORD1/2/3 — written by ChunkSurface.BuildMesh as cornerA/B/C_rgb —
    // as MULTIPLICATIVE tints on top of the height-banded biome colour
    // (`col *= tileTint` at the fragment level). Nearest-corner barycentric
    // selection picks one of the three corner tints per pixel, giving sharp
    // per-cell Voronoi regions instead of soft gradients.
    //
    // Why neutral white for Default + the three point-structure kinds:
    //   • Default          → no tint → cell paints in its natural biome
    //                         colour (grass on land, ocean below sea level).
    //   • Market, Lighthouse, Dock → don't carve, so the terrain wedge IS
    //                         what the player sees around the structure;
    //                         tinting it would also tint the visible ground
    //                         around the building, which reads poorly.
    //   • Road, Building, Field → DO carve in ChunkSurface.BuildMesh; these
    //                         palette entries only paint at cross-chunk
    //                         seams where dvpi == -1 falls through to the
    //                         sentinel (also neutral) — so the entries are
    //                         effectively decorative today. They're kept as
    //                         hints for future "show this kind even on
    //                         foreign corners" extensions.
    //
    // The sentinel index = TileKindPalette.Length - 1 covers `dvpi == -1`
    // foreign corners (the neighbour chunk isn't loaded so we don't know its
    // kind). Neutral so a wedge with a foreign corner doesn't flash a
    // wrong-kind colour at the seam.
    public static class TileKindPalette
    {
        // Order matches TileKind enum byte values 0..6, then the sentinel.
        // Length = 8 = 7 kinds + 1 sentinel. Adding a TileKind requires
        // adding one entry HERE in the same enum order — the build will
        // silently use the sentinel otherwise (which falls back to neutral,
        // so the visual is "no tint" rather than an exception).
        static readonly Color32[] Entries =
        {
            new Color32(255, 255, 255, 255), // 0 Default     — neutral
            new Color32(160, 150, 135, 255), // 1 Road        — warm grey
            new Color32(210, 120, 110, 255), // 2 Building    — soft red
            new Color32(160, 200, 110, 255), // 3 Field       — leafy green
            new Color32(255, 255, 255, 255), // 4 Market      — neutral (structure on top)
            new Color32(255, 255, 255, 255), // 5 Lighthouse  — neutral (structure on top)
            new Color32(255, 255, 255, 255), // 6 Dock        — neutral (structure on top)
            new Color32(255, 255, 255, 255), // 7 sentinel    — foreign corner
        };
        static readonly int SentinelIdx = Entries.Length - 1;

        // Resolve the palette entry for a wedge corner. `primalIdx` is the
        // dual-vertex's primal cell index (= dvpi[cornerCursor + i]) or -1
        // when the corner belongs to a foreign chunk that isn't loaded. `tp`
        // is the chunk's TileProperties array; null means the chunk hasn't
        // been painted (all cells are Default).
        public static Color32 ColorFor(int primalIdx, TileProperty[] tp)
        {
            if (primalIdx < 0 || tp == null || primalIdx >= tp.Length)
                return Entries[SentinelIdx];
            int kind = (int)tp[primalIdx].Kind;
            if ((uint)kind >= (uint)SentinelIdx) return Entries[SentinelIdx];
            return Entries[kind];
        }

        // Pack a Color32 into the Vector3 form the mesh's TEXCOORD1/2/3
        // channels expect. Alpha is dropped because the shader only reads
        // RGB at line 213-215 of TerrainElevation.shader.
        public static Vector3 ColorRgb(Color32 c)
            => new Vector3(c.r / 255f, c.g / 255f, c.b / 255f);
    }
}
