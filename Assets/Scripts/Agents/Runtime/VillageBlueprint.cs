using System.Collections.Generic;
using Agents.Actions;
using Agents.Core;
using Agents.Economy;
using Agents.Goals;
using Agents.Pathfinding;
using Agents.Sensors;
using TerrainGrid;

namespace Agents.Runtime
{
    // One-shot factory for the v1 village blueprint. Constructs the
    // sensor / action / goal lists once at scene start and hands them
    // to every spawned Agent. Same shape as a CrashKonijn AgentBehaviour
    // — just code-driven rather than ScriptableObject-driven.
    //
    // Adding a new role-bearing tile to the demo:
    //   1. Add the TileKind + materials + mesh-builder (Phase A pattern).
    //   2. Add the BuildingRole entry + slot count (Phase B pattern).
    //   3. Add the WorldKey / TargetKey constants (Phase D pattern).
    //   4. Add the new Action + Sensor + (optional) Goal here.
    //   5. Extend VillageSpawner to assign agents into it.
    public static class VillageBlueprint
    {
        public static GoapBlueprint Build(ChunkManager mgr, TilePathCache cache)
        {
            List<GoapSensor> sensors = new List<GoapSensor>
            {
                new TimeOfDaySensor(),
                new NeedsSensor(),
                new InventorySensor(),
                new MoneySensor(),
                new LocationSensor(),
            };

            // Four MoveToTile instances — one per At* key. Each clears
            // the other three At* flags as part of its symbolic effect,
            // so the planner correctly models that an agent can only be
            // at one named location at a time.
            string[] othersField  = { WorldKey.AtBakery, WorldKey.AtMarket, WorldKey.AtHome };
            string[] othersBakery = { WorldKey.AtField,  WorldKey.AtMarket, WorldKey.AtHome };
            string[] othersMarket = { WorldKey.AtField,  WorldKey.AtBakery, WorldKey.AtHome };
            string[] othersHome   = { WorldKey.AtField,  WorldKey.AtBakery, WorldKey.AtMarket };

            List<GoapAction> actions = new List<GoapAction>
            {
                new MoveToTileAction("MoveToField",  WorldKey.AtField,  TargetKey.AssignedField,  othersField,  mgr, cache),
                new MoveToTileAction("MoveToBakery", WorldKey.AtBakery, TargetKey.AssignedBakery, othersBakery, mgr, cache),
                new MoveToTileAction("MoveToMarket", WorldKey.AtMarket, TargetKey.NearestMarket,  othersMarket, mgr, cache),
                new MoveToTileAction("MoveToHome",   WorldKey.AtHome,   TargetKey.Home,           othersHome,   mgr, cache),
                new HarvestWheatAction(),
                new BakeBreadAction(),
                new SellAtMarketAction(Commodity.Wheat),
                new SellAtMarketAction(Commodity.Bread),
                new BuyAtMarketAction(Commodity.Wheat),
                new BuyAtMarketAction(Commodity.Bread),
                new EatAction(),
                new SleepAction(),
            };

            List<GoapGoal> goals = new List<GoapGoal>
            {
                new SurviveGoal(),
                new EarnGoal(),
                new RestockGoal(),
            };

            return new GoapBlueprint(sensors, actions, goals);
        }
    }
}
