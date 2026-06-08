using System;
using System.Collections.Generic;

namespace Agents.Core
{
    // A* over symbolic WorldStates. Pure C# — no Unity, no Agent runtime
    // beyond the IGoapAgent interface for Cost / IsValid dispatch.
    //
    // Nodes: each is a WorldState fingerprint + parent index + the action
    // that produced it from the parent. Reconstruction walks the parent
    // chain from the goal back to the start, emitting actions in order.
    //
    // Heuristic: UnsatisfiedCount(state, target) — the number of target
    // (key, value) pairs not yet present. With v1 actions that mostly
    // toggle a single key per effect, this is close to admissible; an
    // action whose effects satisfy multiple target keys in one shot will
    // make the heuristic overestimate. The cost is loss of *optimality*,
    // not correctness — the planner still finds a valid plan if one
    // exists within MaxNodes.
    //
    // Node cap: 200. v1 plans are 3-7 actions long over a 7-action
    // library, so the actual frontier rarely exceeds 30. The cap is the
    // pathological-input guard. Returns null if exhausted; the caller
    // (Agent.Tick) treats null as "stay on current plan, or drop into
    // Idle if there is none."
    public static class GoapPlanner
    {
        public const int MaxNodes = 200;

        // Plan from `start` (a snapshot of agent.State at plan time) toward
        // the goal. Returns the action sequence as a Queue so the agent
        // can Dequeue() one step at a time; returns an empty Queue if the
        // start already satisfies the goal; returns null if no plan was
        // found.
        public static Queue<GoapAction> Plan(WorldState start, GoapGoal goal,
                                             IReadOnlyList<GoapAction> actions,
                                             IGoapAgent agent)
        {
            if (start == null || goal == null || actions == null) return null;
            IReadOnlyList<KeyValuePair<string, object>> target = MaterializeTarget(goal.TargetState);
            if (start.Satisfies(target)) return new Queue<GoapAction>();

            List<NodeRecord> nodes = new List<NodeRecord>();
            nodes.Add(new NodeRecord { State = start, ParentIdx = -1, ActionTaken = null, G = 0f });

            HashSet<WorldState> closed = new HashSet<WorldState>();
            MinHeap<HeapEntry> heap = new MinHeap<HeapEntry>((a, b) => a.F.CompareTo(b.F));
            heap.Push(new HeapEntry { NodeIdx = 0, F = start.UnsatisfiedCount(target) });

            int explored = 0;
            while (heap.Count > 0 && explored < MaxNodes)
            {
                HeapEntry top = heap.PopMin();
                NodeRecord node = nodes[top.NodeIdx];

                // Stale heap entry (closed-set add fails) ⇒ skip without
                // counting this as work.
                if (!closed.Add(node.State)) continue;
                explored++;

                if (node.State.Satisfies(target))
                    return ReconstructPlan(nodes, top.NodeIdx);

                for (int i = 0; i < actions.Count; i++)
                {
                    GoapAction action = actions[i];
                    if (!node.State.Satisfies(action.Preconditions)) continue;
                    if (!action.IsValid(agent)) continue;

                    WorldState ns = node.State.Clone();
                    ns.Apply(action.Effects);
                    // Cheap pre-check before allocating the heap entry —
                    // the closed set still owns final dedup, so this is
                    // just an early bail-out for the common case.
                    if (closed.Contains(ns)) continue;

                    float newG = node.G + action.Cost(agent);
                    int newH = ns.UnsatisfiedCount(target);

                    int newIdx = nodes.Count;
                    nodes.Add(new NodeRecord
                    {
                        State = ns,
                        ParentIdx = top.NodeIdx,
                        ActionTaken = action,
                        G = newG,
                    });
                    heap.Push(new HeapEntry { NodeIdx = newIdx, F = newG + newH });
                }
            }
            return null;
        }

        // Materialize the target enumerable once so we don't re-iterate
        // (and potentially re-allocate) inside Satisfies on every loop.
        static IReadOnlyList<KeyValuePair<string, object>> MaterializeTarget(
            IEnumerable<KeyValuePair<string, object>> src)
        {
            if (src is IReadOnlyList<KeyValuePair<string, object>> already) return already;
            List<KeyValuePair<string, object>> list = new List<KeyValuePair<string, object>>();
            foreach (KeyValuePair<string, object> kv in src) list.Add(kv);
            return list;
        }

        static Queue<GoapAction> ReconstructPlan(List<NodeRecord> nodes, int endIdx)
        {
            Stack<GoapAction> stack = new Stack<GoapAction>();
            int cur = endIdx;
            // Hard cap so a corrupt parent chain can't infinite-loop.
            int safety = MaxNodes + 8;
            while (cur > 0 && safety-- > 0)
            {
                NodeRecord n = nodes[cur];
                if (n.ActionTaken != null) stack.Push(n.ActionTaken);
                cur = n.ParentIdx;
            }
            return new Queue<GoapAction>(stack);
        }

        struct NodeRecord
        {
            public WorldState State;
            public int ParentIdx;
            public GoapAction ActionTaken;
            public float G;
        }

        struct HeapEntry
        {
            public int NodeIdx;
            public float F;
        }

        // Tiny generic binary min-heap (twin of TilePathfinder's). Inlined
        // here so the planner stays single-file and we don't ship a public
        // utility that demands its own contract.
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
