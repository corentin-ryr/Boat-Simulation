using System;
using System.Collections.Generic;
using System.Text;

namespace Agents.Core
{
    // Symbolic state bag the planner reasons about. Pure data — no Unity
    // types in the public surface so the planner is unit-testable from
    // plain C# fixtures.
    //
    // Modelling choice: STRIPS-style boolean / scalar value comparison
    // rather than numeric reasoning. Preconditions and effects are
    // (key, value) pairs; satisfaction is per-key equality. Numeric
    // thresholds ("Money >= 60") project through a sensor that flips a
    // symbolic bool (`MoneyEnough = true`) when the threshold is crossed.
    //
    // Why a dictionary instead of typed fields:
    //   • New keys cost zero schema change.
    //   • Sensors / actions evolve faster than the planner.
    //   • Matches the CrashKonijn vocabulary — migration to that library
    //     becomes a regex from `Set(WorldKey.X, …)` to `worldData.SetState
    //     <X>(…)`.
    public sealed class WorldState : IEquatable<WorldState>
    {
        // Backing store. Public reads via Get / Has; mutations through Set.
        // Capacity hint matches the v1 key count (~12) so we avoid the
        // initial growth on the first writes.
        readonly Dictionary<string, object> values = new Dictionary<string, object>(16);

        public int Count => values.Count;
        public IEnumerable<KeyValuePair<string, object>> Pairs => values;

        public bool Has(string key) => values.ContainsKey(key);

        public object Get(string key) => values.TryGetValue(key, out object v) ? v : null;

        // Typed lookup. Default(T) is returned if the key is missing OR
        // the stored value isn't the requested type — callers should
        // either guard with Has(key) or trust their sensor contract.
        public T Get<T>(string key)
        {
            if (values.TryGetValue(key, out object v) && v is T t) return t;
            return default;
        }

        public void Set(string key, object value) => values[key] = value;

        public void Remove(string key) => values.Remove(key);

        public void Clear() => values.Clear();

        // Shallow clone — fine for the value types we use (bool, int,
        // string, ValueTuples). If a planner later stores mutable
        // references in the bag, switch to deep clone.
        public WorldState Clone()
        {
            WorldState c = new WorldState();
            foreach (KeyValuePair<string, object> kv in values) c.values[kv.Key] = kv.Value;
            return c;
        }

        // Per-pair containment: does this state hold every (k, v) in
        // `target`? Used by the planner for goal satisfaction and
        // precondition checks. Comparison uses Equals so ValueTuples,
        // strings, ints and bools all behave correctly.
        public bool Satisfies(IEnumerable<KeyValuePair<string, object>> target)
        {
            foreach (KeyValuePair<string, object> kv in target)
            {
                if (!values.TryGetValue(kv.Key, out object v)) return false;
                if (!ValueEquals(v, kv.Value)) return false;
            }
            return true;
        }

        // Overwrite/add the listed pairs. Used by the planner to apply an
        // action's symbolic effects when expanding the search frontier,
        // and by Agent at execution time when an action's Complete() runs.
        public void Apply(IEnumerable<KeyValuePair<string, object>> effects)
        {
            foreach (KeyValuePair<string, object> kv in effects) values[kv.Key] = kv.Value;
        }

        // Count of `target` pairs NOT satisfied by this state. Used as
        // the planner's admissible heuristic: each action can satisfy at
        // most one symbolic precondition at a time (in our v1 schema),
        // so the heuristic never overestimates the remaining plan
        // length.
        public int UnsatisfiedCount(IEnumerable<KeyValuePair<string, object>> target)
        {
            int n = 0;
            foreach (KeyValuePair<string, object> kv in target)
            {
                if (!values.TryGetValue(kv.Key, out object v) || !ValueEquals(v, kv.Value)) n++;
            }
            return n;
        }

        // -------- Equality / hashing for the planner's closed set --------

        public bool Equals(WorldState other)
        {
            if (other == null) return false;
            if (other.values.Count != values.Count) return false;
            foreach (KeyValuePair<string, object> kv in values)
            {
                if (!other.values.TryGetValue(kv.Key, out object v)) return false;
                if (!ValueEquals(v, kv.Value)) return false;
            }
            return true;
        }

        public override bool Equals(object obj) => obj is WorldState ws && Equals(ws);

        // Commutative hash so iteration order doesn't matter. Each pair
        // contributes a salted hash that we XOR into the running total —
        // collisions are tolerable because the closed-set check falls
        // back to Equals when hashes match.
        public override int GetHashCode()
        {
            int h = 0;
            foreach (KeyValuePair<string, object> kv in values)
            {
                int pair = unchecked(kv.Key.GetHashCode() * 397) ^ (kv.Value?.GetHashCode() ?? 0);
                h ^= pair;
            }
            return h;
        }

        static bool ValueEquals(object a, object b)
        {
            if (a == null) return b == null;
            return a.Equals(b);
        }

        public override string ToString()
        {
            if (values.Count == 0) return "(empty)";
            StringBuilder sb = new StringBuilder();
            bool first = true;
            foreach (KeyValuePair<string, object> kv in values)
            {
                if (!first) sb.Append(", ");
                sb.Append(kv.Key).Append('=').Append(kv.Value ?? "null");
                first = false;
            }
            return sb.ToString();
        }
    }
}
