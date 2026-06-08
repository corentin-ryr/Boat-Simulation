namespace Agents.Core
{
    // Marker base for per-agent per-action scratch state. The "actions
    // are stateless, data lives in ActionData" rule borrowed from the
    // CrashKonijn vocabulary — concrete actions never carry mutable
    // fields, so a single action instance can be shared across all
    // agents without bleed-through.
    //
    // Concrete subtypes (Phase F):
    //   • MoveToTileActionData — { TilePath path, int stepIndex }
    //   • HarvestWheatActionData — { float harvestProgressSeconds }
    //   • BakeBreadActionData — { float bakeProgressSeconds }
    //
    // Agents own a Dictionary<Type, ActionData> keyed by the concrete
    // type — IGoapAgent.GetActionData<T>() lazily creates and returns
    // the bucket for the requested type.
    public abstract class ActionData
    {
        // Reset to a fresh state when the action starts. Avoids stale
        // progress timers carrying over from a previous Stop(). Called
        // by Agent before Action.Start().
        public virtual void Reset() { }
    }
}
