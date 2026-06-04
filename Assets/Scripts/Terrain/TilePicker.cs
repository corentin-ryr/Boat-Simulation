using UnityEngine;

namespace TerrainGrid
{
    // Cursor → primal cell resolution, shared by TilePainter (gameplay) and TileSelector
    // (RTS mode). Extracted from the original inline implementation in TilePainter so
    // both pickers walk the same code path: cast a ray from the camera through the
    // cursor, hit the terrain MeshCollider, resolve the hit point to (ChunkCoord, cell
    // index) via the same nearest-centroid lookup the streamer uses, and report the
    // raw hit position for whoever needs to draw a highlight at it.
    //
    // Pure static — no per-instance state. Cheap to call every frame (the only
    // allocations are the local RaycastHit struct and Ray).
    public static class TilePicker
    {
        // Returns true if the ray hits terrain AND the hit chunk has a published primal
        // mirror with centroids. coord/cellIdx identify the picked tile; hitPoint is
        // the raw world-space raycast hit (useful for hover positioning that doesn't
        // need to track the cell centroid exactly).
        public static bool PickUnderCursor(
            Camera cam, ChunkManager mgr, LayerMask mask, float maxDistance,
            out ChunkCoord coord, out int cellIdx, out Vector3 hitPoint)
        {
            coord = default;
            cellIdx = -1;
            hitPoint = default;

            if (cam == null || mgr == null) return false;

            Ray ray = cam.ScreenPointToRay(Input.mousePosition);
            if (!Physics.Raycast(ray, out RaycastHit hit, maxDistance, mask)) return false;
            hitPoint = hit.point;

            // Resolve the hit world position to the chunk it lives in. Reading the hex
            // constants off ChunkManager keeps this aligned with the streamer's coord
            // system without duplicating numbers.
            float hexRadius = mgr.HexRadius;
            int chunkGridSize = mgr.ChunkGridSize;
            coord = ChunkCoord.FromWorldPos(hit.point, hexRadius, chunkGridSize);

            if (!mgr.TryGetPublishedChunk(coord, out PrimalChunk chunk)) return false;
            if (chunk.PrimalCellCentroids == null || chunk.PrimalCellCentroids.Length == 0) return false;

            // Nearest-centroid lookup. Tiles are small (~hexRadius worth of XZ extent)
            // and chunks contain ~50 cells, so this is a few microseconds — no spatial
            // index needed. Y is ignored because elevation can put the centroid above
            // or below the hit point; the painted tile is the one whose XZ footprint
            // owns the cursor.
            int bestIdx = -1;
            float bestSq = float.PositiveInfinity;
            for (int i = 0; i < chunk.PrimalCellCentroids.Length; i++)
            {
                Vector3 c = chunk.PrimalCellCentroids[i];
                float dx = c.x - hit.point.x;
                float dz = c.z - hit.point.z;
                float sq = dx * dx + dz * dz;
                if (sq < bestSq) { bestSq = sq; bestIdx = i; }
            }
            if (bestIdx < 0) return false;

            cellIdx = bestIdx;
            return true;
        }
    }
}
