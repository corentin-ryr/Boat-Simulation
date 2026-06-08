using System.Text;

namespace Agents.Economy
{
    // Per-actor goods storage. Backing store is an int[] indexed by
    // (byte)Commodity — zero per-access allocation, branch-free reads,
    // and no Dictionary boxing. A class rather than a struct so mutations
    // through a reference (e.g. `entry.Stock.Add(...)`) work without
    // chasing-the-copy bugs.
    //
    // Used by:
    //   • Agent — what the agent is carrying.
    //   • BuildingEntry — Field accumulated harvest, Bakery wheat-in /
    //     bread-out, etc.
    //   • MarketStockpile — the market's tradable inventory.
    public class Inventory
    {
        readonly int[] amounts = new int[Commodities.Count];

        public int Get(Commodity c) => amounts[(int)c];

        // Negative deltas are accepted (so callers can pass `Add(c, -1)`
        // when they don't care about the underflow case). For underflow-
        // protected removal use Remove() which returns false on shortage.
        public void Add(Commodity c, int n) => amounts[(int)c] += n;

        // Atomic guard against partial removal: returns false and leaves
        // the inventory untouched when the requested quantity isn't
        // available. Callers should treat this as the contract for any
        // "spend X of commodity Y" intent.
        public bool Remove(Commodity c, int n)
        {
            if (n < 0) return false;
            int idx = (int)c;
            if (amounts[idx] < n) return false;
            amounts[idx] -= n;
            return true;
        }

        public bool Has(Commodity c, int n) => amounts[(int)c] >= n;

        public void Set(Commodity c, int n) => amounts[(int)c] = n < 0 ? 0 : n;

        // Sum across every commodity. Cheap (Commodities.Count is tiny).
        // Used by HUD / debug summaries.
        public int TotalItems()
        {
            int s = 0;
            for (int i = 0; i < amounts.Length; i++) s += amounts[i];
            return s;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            bool first = true;
            for (int i = 0; i < amounts.Length; i++)
            {
                if (amounts[i] <= 0) continue;
                if (!first) sb.Append(", ");
                sb.Append(Commodities.Name((Commodity)i)).Append(':').Append(amounts[i]);
                first = false;
            }
            return first ? "(empty)" : sb.ToString();
        }
    }
}
