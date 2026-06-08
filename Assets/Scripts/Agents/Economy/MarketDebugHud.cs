using UnityEngine;

namespace Agents.Economy
{
    // IMGUI overlay for inspecting + driving Markets while the agent layer
    // isn't yet implemented. Lists every Market the BuildingRegistry knows
    // about and shows Cash / Stock / Price per commodity. Per-Market
    // buttons inject synthetic trades so we can watch the daily price
    // tick respond without waiting for agents.
    //
    // Lives in Economy/ because it's tightly coupled to MarketStockpile +
    // BuildingRegistry. F5 toggles. Disabled by default — the screen real
    // estate gets crowded quickly when many Markets exist.
    public class MarketDebugHud : MonoBehaviour
    {
        public bool drawHud = true;
        public KeyCode toggleKey = KeyCode.F5;

        const float W = 320f;
        const float RowH = 18f;

        void Update()
        {
            if (Input.GetKeyDown(toggleKey)) drawHud = !drawHud;
        }

        void OnGUI()
        {
            if (!drawHud) return;
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;

            var markets = reg.AllOfRole(BuildingRole.Market);

            // Top-of-panel header + a force-day-tick button so we can sim
            // the OnNewDay price adjustment without sleeping a full in-game
            // day.
            float x = 10f;
            float y = 70f;
            int rows = 2 + markets.Count * (3 + Commodities.Count);
            float h = (rows + 2) * RowH;
            GUI.Box(new Rect(x, y, W, h), "");
            GUI.Label(new Rect(x + 6f, y + 2f, W - 12f, RowH),
                      $"Markets ({markets.Count})  [{toggleKey} hide]");
            if (GUI.Button(new Rect(x + 6f, y + 2f + RowH, 140f, RowH),
                           "Force daily tick"))
                ForceDailyTickAllMarkets();

            float ry = y + 2f + 2 * RowH;
            for (int i = 0; i < markets.Count; i++)
            {
                BuildingEntry e = markets[i];
                MarketStockpile m = e.Market;
                if (m == null) continue;

                GUI.Label(new Rect(x + 6f, ry, W - 12f, RowH),
                          $"#{i} {e.Coord}:{e.CellIndex}  cash=${m.Cash}");
                ry += RowH;

                foreach (Commodity c in Commodities.All)
                {
                    int stock = m.Stock.Get(c);
                    float price = m.GetPrice(c);
                    GUI.Label(new Rect(x + 16f, ry, 160f, RowH),
                              $"{Commodities.Name(c)}: {stock} @ ${price:0.0}");
                    if (GUI.Button(new Rect(x + 180f, ry, 60f, RowH), "+1 sell"))
                        m.TryBuyFrom(c, 1, out _);
                    if (GUI.Button(new Rect(x + 246f, ry, 60f, RowH), "-1 buy"))
                        m.TrySellTo(c, 1, 999, out _);
                    ry += RowH;
                }
                ry += 4f;
            }
        }

        // Walks every Market and applies one PriceModel.DailyTick. Useful
        // for verifying that paint-a-Field-and-watch-prices-drop loop
        // without waiting for the GameClock to wrap a day.
        void ForceDailyTickAllMarkets()
        {
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null) return;
            foreach (BuildingEntry e in reg.AllOfRole(BuildingRole.Market))
                if (e.Market != null) PriceModel.DailyTick(e.Market);
        }
    }
}
