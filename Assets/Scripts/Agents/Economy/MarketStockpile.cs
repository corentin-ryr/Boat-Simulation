using System.Collections.Generic;
using UnityEngine;

namespace Agents.Economy
{
    // Per-Market trading state. One instance per Market tile, owned by the
    // corresponding BuildingEntry. Markets are buyers AND sellers — a
    // farmer sells wheat to the market (Market.Cash → farmer; wheat →
    // Market.Stock); a citizen buys bread from the market (citizen cash
    // → Market.Cash; bread → citizen).
    //
    // Prices and cash float are per-Market (no global price for v1) so
    // distant markets can drift apart, which is the seed for inter-Market
    // trade routes when we add them.
    public class MarketStockpile
    {
        public Inventory Stock { get; } = new Inventory();
        public int Cash;

        readonly float[] prices = new float[Commodities.Count];

        // Stock at the start of the in-game day. Snapshotted by
        // PriceModel.DailyTick after it adjusts prices, so the next day's
        // delta is measured against the right baseline.
        readonly int[] lastDayStock = new int[Commodities.Count];

        // Rolling history shown by the Market HUD. Bounded by RecentTradeCap
        // so even hundreds of trades a day don't grow without bound.
        const int RecentTradeCap = 32;
        readonly Queue<TradeRecord> recentTrades = new Queue<TradeRecord>();
        public IReadOnlyCollection<TradeRecord> RecentTrades => recentTrades;

        public MarketStockpile()
        {
            for (int i = 0; i < prices.Length; i++)
                prices[i] = PriceModel.InitialPrice((Commodity)i);
            Cash = PriceModel.InitialMarketCash;
            // Snapshot the (currently empty) start-of-day stock so the
            // first OnNewDay tick has a real baseline rather than a NRE.
            SnapshotEndOfDay();
        }

        public float GetPrice(Commodity c) => prices[(int)c];

        public void SetPrice(Commodity c, float p)
            => prices[(int)c] = Mathf.Clamp(p, PriceModel.MinPrice, PriceModel.MaxPrice);

        public int LastDayStockOf(Commodity c) => lastDayStock[(int)c];

        public void SnapshotEndOfDay()
        {
            for (int i = 0; i < lastDayStock.Length; i++)
                lastDayStock[i] = Stock.Get((Commodity)i);
        }

        // Market buys `qty` of `c` from a seller (e.g. a farmer offloading
        // wheat). Returns false (no state change) if the market can't
        // afford it. `agentInventoryCheck` is the caller's promise that
        // the goods exist — Stockpile only tracks its own side of the
        // ledger.
        public bool TryBuyFrom(Commodity c, int qty, out int totalPaid)
        {
            totalPaid = 0;
            if (qty <= 0) return false;
            float pricePer = GetPrice(c);
            int payout = Mathf.FloorToInt(pricePer * qty);
            if (Cash < payout) return false;
            Stock.Add(c, qty);
            Cash -= payout;
            totalPaid = payout;
            Record(new TradeRecord(c, qty, payout, TradeDirection.MarketBought));
            return true;
        }

        // Customer pays `agentCashAvailable` for up to `qty` of `c`. Returns
        // false if the market is out of stock, the agent can't afford even
        // the full quantity, or qty is invalid. Failure leaves both sides
        // untouched.
        public bool TrySellTo(Commodity c, int qty, int agentCashAvailable, out int totalCost)
        {
            totalCost = 0;
            if (qty <= 0) return false;
            if (!Stock.Has(c, qty)) return false;
            float pricePer = GetPrice(c);
            int cost = Mathf.CeilToInt(pricePer * qty);
            if (agentCashAvailable < cost) return false;
            Stock.Remove(c, qty);
            Cash += cost;
            totalCost = cost;
            Record(new TradeRecord(c, qty, cost, TradeDirection.MarketSold));
            return true;
        }

        void Record(TradeRecord t)
        {
            recentTrades.Enqueue(t);
            while (recentTrades.Count > RecentTradeCap) recentTrades.Dequeue();
        }
    }

    public enum TradeDirection : byte
    {
        // Market is the buyer (agent sold to market). Reads as cash leaving
        // the market and goods entering its stockpile.
        MarketBought = 0,
        // Market is the seller (agent bought from market). Cash in, goods out.
        MarketSold   = 1,
    }

    public readonly struct TradeRecord
    {
        public readonly Commodity Commodity;
        public readonly int Quantity;
        public readonly int Money;
        public readonly TradeDirection Direction;
        public TradeRecord(Commodity c, int qty, int money, TradeDirection dir)
        {
            Commodity = c; Quantity = qty; Money = money; Direction = dir;
        }
    }
}
