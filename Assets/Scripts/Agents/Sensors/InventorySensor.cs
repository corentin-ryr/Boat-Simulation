using Agents.Core;
using Agents.Economy;

namespace Agents.Sensors
{
    // Projects Inventory commodity counts into the bool WorldKeys the
    // planner uses. Adding a tracked commodity = one new line here.
    public sealed class InventorySensor : GoapSensor
    {
        public override string Name => "Inventory";
        public override SenseMode Mode => SenseMode.Always;

        public override void Sense(IGoapAgent agent)
        {
            agent.State.Set(WorldKey.HasWheat, agent.Inventory.Has(Commodity.Wheat, 1));
            agent.State.Set(WorldKey.HasBread, agent.Inventory.Has(Commodity.Bread, 1));
        }
    }
}
