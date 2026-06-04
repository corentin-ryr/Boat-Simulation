using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid
{
    // The full registry of TileDefinitions for one scene. Inspector-wired onto
    // ChunkManager (one catalog per ChunkManager). Read in two patterns:
    //   - Sequential iteration (TilePalette HUD, ChunkSurface EnsurePresence)
    //   - Kind → Definition lookup (TilePlacement IsValid, ChunkSurface BuildMesh
    //     carve gate)
    //
    // The lookup table is rebuilt on demand the first time `Get` is called (and
    // any time `definitions` changes — detected via length mismatch). Cheap
    // enough that we don't need an explicit invalidation API.
    [CreateAssetMenu(menuName = "Terrain/Tile Catalog", fileName = "TileCatalog")]
    public class TileCatalog : ScriptableObject
    {
        [Tooltip("All tile kinds available to the player in this scene. Order is the " +
                 "palette HUD order and the cycle-key order. The Default entry should " +
                 "come first so an unpainted scene reads naturally in the palette.")]
        public List<TileDefinition> definitions = new List<TileDefinition>();

        // Built lazily from `definitions`. Reset to null when the source list
        // size changes (catches Inspector edits between calls).
        Dictionary<TileKind, TileDefinition> byKind;
        int lastBuiltCount = -1;

        public IReadOnlyList<TileDefinition> Definitions => definitions;

        // Returns the definition for this kind, or null if the kind has no
        // entry. Most callers should handle null gracefully (e.g. by treating
        // it as Default + no presenter).
        public TileDefinition Get(TileKind kind)
        {
            EnsureIndex();
            return byKind.TryGetValue(kind, out TileDefinition def) ? def : null;
        }

        // True when the terrain mesh should skip cells of this kind's wedge in
        // favour of the kind's presenter. Centralised here so BuildMesh and the
        // presenter don't disagree about which kinds carve.
        public bool IsCarvingKind(TileKind kind)
        {
            TileDefinition def = Get(kind);
            return def != null && def.carvesTerrain;
        }

        // Resolve the Material a presenter should use for the given kind. Null
        // if the kind has no definition or no material wired — caller decides
        // what to do (Unity falls back to magenta, which is loud enough to
        // notice during scene setup).
        public Material MaterialFor(TileKind kind)
        {
            TileDefinition def = Get(kind);
            return def?.material;
        }

        void EnsureIndex()
        {
            if (byKind != null && lastBuiltCount == definitions.Count) return;
            byKind = new Dictionary<TileKind, TileDefinition>(definitions.Count);
            foreach (TileDefinition def in definitions)
            {
                if (def == null) continue;
                // Last entry wins on duplicate Kind — log so a misconfigured
                // scene is noisy at startup rather than silently picking one.
                if (byKind.ContainsKey(def.kind))
                    Debug.LogWarning($"TileCatalog '{name}': duplicate definition for kind {def.kind} — using '{def.name}'.");
                byKind[def.kind] = def;
            }
            lastBuiltCount = definitions.Count;
        }
    }
}
