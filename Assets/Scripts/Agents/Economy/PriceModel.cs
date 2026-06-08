namespace Agents.Economy
{
    // Per-Market price tuning. Held as a static helper rather than a
    // ScriptableObject because it's a small numeric table — code-driven
    // matches the rest of the agent layer.
    //
    // Daily price dynamics: every Market's prices nudge up when stock
    // went DOWN over the day (undersupplied → harder to buy) and down
    // when stock went UP (oversupplied → easier to clear). Clamped to
    // [MinPrice, MaxPrice] so a runaway dynamic can't print arbitrarily
    // large numbers. Not a real auction — just a heuristic that produces
    // the headline behaviour ("paint a second Field → bread cheaper").
    public static class PriceModel
    {
        public const float MinPrice = 1f;
        public const float MaxPrice = 200f;

        // Multiplicative adjustment per day. Smaller than 10 % so prices
        // drift slowly enough to be readable. Asymmetric (up faster than
        // down) so undersupply hurts visibly before the market self-
        // corrects.
        public const float AdjustUp   = 1.06f;
        public const float AdjustDown = 0.95f;

        public static float InitialPrice(Commodity c) => c switch
        {
            Commodity.Wheat => 5f,
            Commodity.Bread => 12f,
            _               => 10f,
        };

        // Initial cash float seeded into every new Market. High enough
        // that the Market can buy the first few rounds of harvest before
        // any customers spend; low enough that runaway accumulation is
        // self-limiting.
        public const int InitialMarketCash = 200;

        // Move each commodity's price according to its overnight stock
        // delta, then snapshot the new "start-of-day" stock so the next
        // tick has a baseline. Called by BuildingRegistry on
        // GameClock.OnNewDay.
        public static void DailyTick(MarketStockpile m)
        {
            foreach (Commodity c in Commodities.All)
            {
                int now = m.Stock.Get(c);
                int prev = m.LastDayStockOf(c);
                float price = m.GetPrice(c);
                int delta = now - prev;
                if (delta > 0) price *= AdjustDown;     // oversupply
                else if (delta < 0) price *= AdjustUp;  // undersupply
                m.SetPrice(c, price);
            }
            m.SnapshotEndOfDay();
        }
    }
}
