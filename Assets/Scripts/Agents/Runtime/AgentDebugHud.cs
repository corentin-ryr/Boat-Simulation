using System.Text;
using UnityEngine;

namespace Agents.Runtime
{
    // IMGUI overlay listing every live Agent + their internal state.
    // F3 toggles. The primary verification tool for Phase E and F — the
    // GOAP loop is opaque from the outside, so the only way to know
    // what an agent is thinking is to expose it here.
    //
    // Columns: name, goal, current action, plan length, money, hunger,
    // energy, inventory summary, at-tile.
    public class AgentDebugHud : MonoBehaviour
    {
        public bool drawHud = true;
        public KeyCode toggleKey = KeyCode.F3;

        const float W = 720f;
        const float RowH = 18f;
        const float HeaderH = 20f;

        // Scroll position so the panel doesn't overflow when there are
        // many agents in the scene.
        Vector2 scroll;

        readonly StringBuilder sb = new StringBuilder(128);

        void Update()
        {
            if (Input.GetKeyDown(toggleKey)) drawHud = !drawHud;
        }

        void OnGUI()
        {
            if (!drawHud) return;
            AgentRegistry reg = AgentRegistry.Instance;
            if (reg == null) return;

            float x = 10f;
            float y = 130f;
            float h = HeaderH + RowH + Mathf.Min(reg.Count, 30) * RowH + 8f;
            GUI.Box(new Rect(x, y, W, h), "");

            GUI.Label(new Rect(x + 6f, y + 2f, W - 12f, HeaderH),
                      $"Agents ({reg.Count})  [{toggleKey} hide]");

            // Column header.
            float rowY = y + 2f + HeaderH;
            GUILayout.BeginArea(new Rect(x + 6f, rowY, W - 12f, h - HeaderH));
            DrawHeaderRow();
            scroll = GUILayout.BeginScrollView(scroll, GUILayout.Width(W - 16f),
                                               GUILayout.Height(h - HeaderH - RowH - 8f));
            for (int i = 0; i < reg.All.Count; i++)
            {
                Agent a = reg.All[i];
                if (a == null) continue;
                DrawAgentRow(a);
            }
            GUILayout.EndScrollView();
            GUILayout.EndArea();
        }

        void DrawHeaderRow()
        {
            GUILayout.BeginHorizontal();
            GUILayout.Label("name",       GUILayout.Width(80));
            GUILayout.Label("goal",       GUILayout.Width(90));
            GUILayout.Label("action",     GUILayout.Width(120));
            GUILayout.Label("plan",       GUILayout.Width(40));
            GUILayout.Label("$",          GUILayout.Width(40));
            GUILayout.Label("hunger",     GUILayout.Width(50));
            GUILayout.Label("energy",     GUILayout.Width(50));
            GUILayout.Label("inv",        GUILayout.Width(120));
            GUILayout.Label("at",         GUILayout.Width(100));
            GUILayout.EndHorizontal();
        }

        void DrawAgentRow(Agent a)
        {
            GUILayout.BeginHorizontal();
            GUILayout.Label(a.DisplayName, GUILayout.Width(80));
            GUILayout.Label(a.CurrentGoal?.Name ?? "—", GUILayout.Width(90));
            GUILayout.Label(a.CurrentAction?.Name ?? "—", GUILayout.Width(120));
            GUILayout.Label((a.CurrentPlan?.Count ?? 0).ToString(), GUILayout.Width(40));
            GUILayout.Label(a.Money.ToString(), GUILayout.Width(40));
            GUILayout.Label($"{a.Needs.Hunger:0}", GUILayout.Width(50));
            GUILayout.Label($"{a.Needs.Energy:0}", GUILayout.Width(50));
            GUILayout.Label(a.Inventory.ToString(), GUILayout.Width(120));
            GUILayout.Label(a.AtTile.HasValue
                ? $"{a.AtTile.Value.Coord}:{a.AtTile.Value.CellIndex}"
                : "—", GUILayout.Width(100));
            GUILayout.EndHorizontal();
        }
    }
}
