using System;
using System.Collections.Generic;
using TerrainGrid;
using UnityEngine;

namespace Agents.Pathfinding
{
    // A* over the primal cell graph, cross-chunk aware via PrimalChunk's
    // CSR neighbour tables (PrimalEdgeOffsets + PrimalNeighborChunkDQ/DR +
    // PrimalNeighborPolyIdx). Same traversal pattern TilePlacement uses
    // for adjacency tests, just walked greedily for shortest-path search.
    //
    // Walkability model:
    //   • Cells below sea level are blocked unless the kind is Dock or
    //     Lighthouse (both ride on the harbour platform at Y = 2).
    //   • Roofed building kinds (Building / House / Bakery) are blocked
    //     except when the cell IS the destination — agents can step "in
    //     through the door" but never route THROUGH a neighbour's house.
    //   • Everything else is walkable.
    //
    // Cost model:
    //   step_cost = 0.5 × (cellCost[from] + cellCost[to]) × euclideanDist
    //   cellCost: Road = 0.5, Field = 1.2, others = 1.0.
    //   Heuristic = MinCellCost × distance-to-target (admissible &
    //   consistent because MinCellCost is the smallest multiplier any
    //   cell can have).
    //
    // The planner is bounded by MaxNodesExplored so an unreachable target
    // fails in microseconds instead of stalling the frame. Failure returns
    // null; the caller (an Agent's plan check) treats null as "drop this
    // plan, replan or fall back."
    public static class TilePathfinder
    {
        public const int MaxNodesExplored = 4000;
        public const float MinCellCost = 0.5f;
        const float SeaLevelY = 0f;
        const float SeaLevelEpsilon = 0.001f;

        public static TilePath FindPath(ChunkCoord fromCoord, int fromCell,
                                        ChunkCoord toCoord, int toCell,
                                        ChunkManager mgr)
        {
            if (mgr == null) return null;
            if (!mgr.TryGetPublishedChunk(fromCoord, out PrimalChunk fromChunk)) return null;
            if (!mgr.TryGetPublishedChunk(toCoord, out PrimalChunk toChunk)) return null;
            Vector3[] fromCentroids = fromChunk.PrimalCellCentroids;
            Vector3[] toCentroids = toChunk.PrimalCellCentroids;
            if (fromCentroids == null || fromCell < 0 || fromCell >= fromCentroids.Length) return null;
            if (toCentroids == null || toCell < 0 || toCell >= toCentroids.Length) return null;

            // Same-cell trivial path.
            if (fromCoord.Equals(toCoord) && fromCell == toCell)
            {
                List<TileStep> trivial = new List<TileStep>
                {
                    new TileStep(fromCoord, fromCell, fromCentroids[fromCell]),
                };
                return new TilePath(trivial, 0f);
            }

            // If the target cell itself isn't walkable as a destination,
            // bail. (E.g. a Lighthouse cell can host a tower — still
            // walkable for routing because the player can stand on the
            // platform; but a Building cell only matters if the agent
            // wants to BE there.)
            TileKind toKind = KindAt(toChunk, toCell);
            if (!IsWalkable(toKind, toCentroids[toCell].y, isTarget: true)) return null;

            Vector3 targetPos = toCentroids[toCell];

            // A* working state. Tuples as composite keys — heap nodes carry
            // their own coord/cell so we don't allocate a per-node object.
            Dictionary<(ChunkCoord, int), float> bestG = new Dictionary<(ChunkCoord, int), float>();
            Dictionary<(ChunkCoord, int), (ChunkCoord, int)> cameFrom
                = new Dictionary<(ChunkCoord, int), (ChunkCoord, int)>();
            HashSet<(ChunkCoord, int)> closed = new HashSet<(ChunkCoord, int)>();
            MinHeap<OpenNode> heap = new MinHeap<OpenNode>((a, b) => a.F.CompareTo(b.F));

            (ChunkCoord, int) startKey = (fromCoord, fromCell);
            bestG[startKey] = 0f;
            float startH = MinCellCost * (fromCentroids[fromCell] - targetPos).magnitude;
            heap.Push(new OpenNode { Coord = fromCoord, Cell = fromCell, G = 0f, F = startH });

            int explored = 0;
            while (heap.Count > 0 && explored < MaxNodesExplored)
            {
                OpenNode cur = heap.PopMin();
                (ChunkCoord, int) curKey = (cur.Coord, cur.Cell);
                if (!closed.Add(curKey)) continue;  // stale heap entry (a better G was found)
                explored++;

                if (cur.Coord.Equals(toCoord) && cur.Cell == toCell)
                    return ReconstructPath(cameFrom, mgr, startKey, curKey, cur.G);

                if (!mgr.TryGetPublishedChunk(cur.Coord, out PrimalChunk curChunk)) continue;
                int[] off = curChunk.PrimalEdgeOffsets;
                int[] dq = curChunk.PrimalNeighborChunkDQ;
                int[] dr = curChunk.PrimalNeighborChunkDR;
                int[] nbIdx = curChunk.PrimalNeighborPolyIdx;
                if (off == null || dq == null || dr == null || nbIdx == null) continue;
                if (cur.Cell + 1 >= off.Length) continue;

                Vector3[] curCentroids = curChunk.PrimalCellCentroids;
                if (curCentroids == null || cur.Cell >= curCentroids.Length) continue;
                Vector3 curPos = curCentroids[cur.Cell];
                TileKind curKind = KindAt(curChunk, cur.Cell);
                float curCellCost = CellCost(curKind);

                int start = off[cur.Cell];
                int end = off[cur.Cell + 1];
                for (int j = start; j < end; j++)
                {
                    int nbCell = nbIdx[j];
                    if (nbCell < 0) continue;
                    ChunkCoord nbCoord = new ChunkCoord(cur.Coord.q + dq[j], cur.Coord.r + dr[j]);
                    PrimalChunk nbChunk;
                    if (dq[j] == 0 && dr[j] == 0) nbChunk = curChunk;
                    else if (!mgr.TryGetPublishedChunk(nbCoord, out nbChunk)) continue;

                    Vector3[] nbCentroids = nbChunk.PrimalCellCentroids;
                    if (nbCentroids == null || nbCell >= nbCentroids.Length) continue;

                    bool isTargetCell = nbCoord.Equals(toCoord) && nbCell == toCell;
                    TileKind nbKind = KindAt(nbChunk, nbCell);
                    Vector3 nbPos = nbCentroids[nbCell];
                    if (!IsWalkable(nbKind, nbPos.y, isTargetCell)) continue;

                    (ChunkCoord, int) nbKey = (nbCoord, nbCell);
                    if (closed.Contains(nbKey)) continue;

                    float nbCellCost = CellCost(nbKind);
                    float stepDist = (nbPos - curPos).magnitude;
                    float stepCost = 0.5f * (curCellCost + nbCellCost) * stepDist;
                    float tentativeG = cur.G + stepCost;

                    if (bestG.TryGetValue(nbKey, out float existing) && existing <= tentativeG) continue;
                    bestG[nbKey] = tentativeG;
                    cameFrom[nbKey] = curKey;
                    float h = MinCellCost * (nbPos - targetPos).magnitude;
                    heap.Push(new OpenNode { Coord = nbCoord, Cell = nbCell, G = tentativeG, F = tentativeG + h });
                }
            }
            return null;
        }

        // Walk back through the parent map and emit the steps in order.
        // Each step's WorldPos is the centroid as-published right now —
        // the caller can drive the agent through these positions; if a
        // chunk later unloads, the step's position is stale but the
        // (coord, cell) handle stays meaningful.
        static TilePath ReconstructPath(Dictionary<(ChunkCoord, int), (ChunkCoord, int)> cameFrom,
                                        ChunkManager mgr,
                                        (ChunkCoord, int) startKey,
                                        (ChunkCoord, int) endKey,
                                        float totalCost)
        {
            List<TileStep> steps = new List<TileStep>();
            (ChunkCoord, int) cur = endKey;
            // Hard cap on chain walk so a corrupt parent map can't infinite-loop.
            int safety = MaxNodesExplored + 16;
            while (safety-- > 0)
            {
                Vector3 pos = LookupCentroid(mgr, cur.Item1, cur.Item2);
                steps.Add(new TileStep(cur.Item1, cur.Item2, pos));
                if (cur.Item1.Equals(startKey.Item1) && cur.Item2 == startKey.Item2) break;
                if (!cameFrom.TryGetValue(cur, out (ChunkCoord, int) prev)) break;
                cur = prev;
            }
            steps.Reverse();
            return new TilePath(steps, totalCost);
        }

        static Vector3 LookupCentroid(ChunkManager mgr, ChunkCoord coord, int cell)
        {
            if (!mgr.TryGetPublishedChunk(coord, out PrimalChunk pc)) return Vector3.zero;
            if (pc.PrimalCellCentroids == null || cell >= pc.PrimalCellCentroids.Length) return Vector3.zero;
            return pc.PrimalCellCentroids[cell];
        }

        static TileKind KindAt(PrimalChunk chunk, int cell)
        {
            TileProperty[] tp = chunk.TileProperties;
            if (tp == null || cell < 0 || cell >= tp.Length) return TileKind.Default;
            return tp[cell].Kind;
        }

        public static bool IsWalkable(TileKind kind, float centroidY, bool isTarget)
        {
            // Below sea level: only platforms ride above the water.
            if (centroidY < SeaLevelY - SeaLevelEpsilon)
                return kind == TileKind.Dock || kind == TileKind.Lighthouse;

            // Solid roofed structures only enterable as a destination.
            if (kind == TileKind.Building || kind == TileKind.House || kind == TileKind.Bakery)
                return isTarget;

            return true;
        }

        public static float CellCost(TileKind kind) => kind switch
        {
            TileKind.Road  => 0.5f,
            TileKind.Field => 1.2f,
            _              => 1.0f,
        };

        struct OpenNode
        {
            public ChunkCoord Coord;
            public int Cell;
            public float G;
            public float F;
        }

        // Tiny binary min-heap. Inlined here because it's only used by this
        // pathfinder and we don't want a public utility floating around the
        // codebase that demands its own contract. Generic so we can swap in
        // other node types later without rewriting it.
        sealed class MinHeap<T>
        {
            readonly List<T> items = new List<T>();
            readonly Comparison<T> cmp;
            public MinHeap(Comparison<T> cmp) { this.cmp = cmp; }
            public int Count => items.Count;
            public void Push(T item) { items.Add(item); SiftUp(items.Count - 1); }
            public T PopMin()
            {
                T top = items[0];
                int last = items.Count - 1;
                items[0] = items[last];
                items.RemoveAt(last);
                if (items.Count > 0) SiftDown(0);
                return top;
            }
            void SiftUp(int i)
            {
                while (i > 0)
                {
                    int p = (i - 1) / 2;
                    if (cmp(items[i], items[p]) < 0)
                    {
                        (items[i], items[p]) = (items[p], items[i]);
                        i = p;
                    }
                    else break;
                }
            }
            void SiftDown(int i)
            {
                int n = items.Count;
                while (true)
                {
                    int l = 2 * i + 1;
                    int r = 2 * i + 2;
                    int m = i;
                    if (l < n && cmp(items[l], items[m]) < 0) m = l;
                    if (r < n && cmp(items[r], items[m]) < 0) m = r;
                    if (m == i) break;
                    (items[m], items[i]) = (items[i], items[m]);
                    i = m;
                }
            }
        }
    }
}
