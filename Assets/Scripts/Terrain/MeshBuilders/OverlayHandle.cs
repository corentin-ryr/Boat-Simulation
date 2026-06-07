using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Per-chunk per-kind overlay state. ChunkSurface owns one instance per
    // (chunk, TileKind) pair the chunk has ever had — stored sparsely in
    // Presence.Overlays — and the dispatch helper RebuildOverlay populates /
    // tears down its GameObject + Mesh + MeshCollider through this struct.
    //
    // Plain class so reads/writes don't need ref tracking; one allocation per
    // overlay lifetime is negligible (chunks live for many frames once
    // resident, and the warm-presence cache reuses handles across revivals).
    public class OverlayHandle
    {
        public GameObject Go;
        public MeshFilter Mf;
        public MeshRenderer Mr;
        public MeshCollider Mc;          // null until requiresCollider && wantCollider
        public Mesh Mesh;
        public int BuiltFromVersion = -1;
    }
}
