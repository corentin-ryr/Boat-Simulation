using System.Collections.Generic;

namespace Agents.Pathfinding
{
    // Sequence of primal-cell centroids forming a walkable route. Produced
    // by TilePathfinder.FindPath; consumed by MoveToTileAction (Phase F) and
    // by the path-debug overlay. Immutable after construction so multiple
    // agents could share a cached path safely (TilePathCache keys are
    // (from, to) so this is the common case).
    public sealed class TilePath
    {
        public IReadOnlyList<TileStep> Steps { get; }
        public float TotalCost { get; }
        public int Length => Steps.Count;
        public TileStep this[int i] => Steps[i];

        public TilePath(IReadOnlyList<TileStep> steps, float totalCost)
        {
            Steps = steps;
            TotalCost = totalCost;
        }

        public TileStep Start => Steps[0];
        public TileStep End => Steps[Steps.Count - 1];
    }
}
