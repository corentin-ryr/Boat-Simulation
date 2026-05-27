using System.Collections.Generic;

namespace TerrainGrid
{
    // Persistent store for chunk snapshots. Unloading frees the heavy runtime objects
    // (GameObjects/meshes and the live Vertex/Polygon graph) while the chunk's relaxed state
    // survives here, so a reload restores it exactly. Backends are swappable (memory now,
    // disk later) behind this interface.
    public interface IChunkStore
    {
        bool Has(ChunkCoord coord);
        ChunkSnapshot Load(ChunkCoord coord);
        void Save(ChunkCoord coord, ChunkSnapshot snapshot);
    }

    // Session-lifetime in-memory store. Holds compact snapshots (primitive arrays), which are
    // far lighter than retaining the live chunks, and lets the render layer free its meshes.
    // RAM grows with explored area — swap for a disk-backed store for bounded/cross-session use.
    public class MemoryChunkStore : IChunkStore
    {
        readonly Dictionary<ChunkCoord, ChunkSnapshot> store = new();

        public bool Has(ChunkCoord coord) => store.ContainsKey(coord);
        public ChunkSnapshot Load(ChunkCoord coord) => store[coord];
        public void Save(ChunkCoord coord, ChunkSnapshot snapshot) => store[coord] = snapshot;
    }
}
