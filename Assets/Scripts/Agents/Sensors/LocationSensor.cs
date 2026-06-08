using Agents.Core;
using Agents.Economy;
using TerrainGrid;
using UnityEngine;

namespace Agents.Sensors
{
    // Projects "am I at my assigned <X>?" booleans from agent.Position
    // against the four target-key tuples stored in WorldState. A small
    // proximity threshold accounts for the kinematic tween in
    // MoveToTileAction stopping near (not exactly at) the target
    // centroid.
    public sealed class LocationSensor : GoapSensor
    {
        public override string Name => "Location";
        public override SenseMode Mode => SenseMode.Always;

        // Squared distance threshold below which an agent counts as "at"
        // a tile. 1.5 m fits comfortably inside one dual cell (~2 m
        // across at hexRadius=1).
        const float ArrivalRadius = 1.5f;
        const float ArrivalRadiusSq = ArrivalRadius * ArrivalRadius;

        public override void Sense(IGoapAgent agent)
        {
            SenseAt(agent, WorldKey.AtField,  TargetKey.AssignedField);
            SenseAt(agent, WorldKey.AtBakery, TargetKey.AssignedBakery);
            SenseAt(agent, WorldKey.AtMarket, TargetKey.NearestMarket);
            SenseAt(agent, WorldKey.AtHome,   TargetKey.Home);
        }

        static void SenseAt(IGoapAgent agent, string atKey, string targetKey)
        {
            if (!agent.State.Has(targetKey)) { agent.State.Set(atKey, false); return; }
            var tup = agent.State.Get<(ChunkCoord Coord, int CellIndex)>(targetKey);
            BuildingRegistry reg = BuildingRegistry.Instance;
            if (reg == null || !reg.TryGet(tup.Coord, tup.CellIndex, out BuildingEntry e))
            {
                agent.State.Set(atKey, false);
                return;
            }
            // XZ-only distance: agents stand on the ground, the centroid
            // is at terrain Y. Comparing in 3D would mark "at Bakery" as
            // false the moment the agent's capsule centre is taller than
            // the centroid.
            Vector3 p = agent.Position; p.y = 0f;
            Vector3 t = e.WorldPos;     t.y = 0f;
            agent.State.Set(atKey, (p - t).sqrMagnitude < ArrivalRadiusSq);
        }
    }
}
