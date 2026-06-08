using System.Collections.Generic;
using Agents.Core;

namespace Agents.Runtime
{
    // Shared (sensors, actions, goals) bundle that every Agent runs
    // against. Stateless concrete instances are safe to share — the
    // CrashKonijn discipline that actions own no per-agent fields means
    // one BakeBreadAction object can serve all 100 Bakers concurrently.
    //
    // Constructed once at scene start (VillageSpawner in Phase F), assigned
    // to each Agent. Distinct blueprints can coexist if we ever want, say,
    // a Caravan agent with different goals — for v1 there's one blueprint.
    public sealed class GoapBlueprint
    {
        public IReadOnlyList<GoapSensor> Sensors { get; }
        public IReadOnlyList<GoapAction> Actions { get; }
        public IReadOnlyList<GoapGoal>   Goals   { get; }

        public GoapBlueprint(IReadOnlyList<GoapSensor> sensors,
                             IReadOnlyList<GoapAction> actions,
                             IReadOnlyList<GoapGoal> goals)
        {
            Sensors = sensors ?? System.Array.Empty<GoapSensor>();
            Actions = actions ?? System.Array.Empty<GoapAction>();
            Goals   = goals   ?? System.Array.Empty<GoapGoal>();
        }
    }
}
