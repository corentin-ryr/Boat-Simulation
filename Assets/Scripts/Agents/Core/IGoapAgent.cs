using Agents.Economy;
using TerrainGrid;
using UnityEngine;

namespace Agents.Core
{
    // The contract the planner + actions use to talk to an agent. The
    // Phase E Agent MonoBehaviour implements this; planner-side unit
    // tests can supply a fake IGoapAgent that returns canned values.
    //
    // Read-only properties from the planner's perspective: State,
    // Position, Home, Job, Inventory. Mutable for Money so trading
    // actions can settle payouts. Inventory mutations route through
    // its own typed API.
    //
    // ActionData lives on the agent, retrieved by type. The planner
    // doesn't touch it — only concrete actions do, during Start/Perform/
    // Complete. Keeps the planner uninterested in execution scratch.
    public interface IGoapAgent
    {
        // The agent's symbolic state — sensors write here, planner reads.
        WorldState State { get; }

        // World-space position used by sensors (At-tile projection,
        // nearest-Market lookup) and by movement actions for their
        // tween. Mutable so kinematic actions can drive it directly;
        // a future physics-driven mover would expose a different
        // wrapper.
        Vector3 Position { get; set; }

        // Inventory of carryable Commodities. Mutated atomically via
        // the Inventory class's own methods.
        Inventory Inventory { get; }

        // Cash on hand. Settled by trading actions.
        int Money { get; set; }

        // Cells the agent treats as targets when planning. The tuple is
        // (chunk coord, primal cell index). Nullable so agents that
        // haven't been assigned a role yet aren't forced into bogus
        // home/job values.
        (ChunkCoord Coord, int CellIndex)? AtTile { get; set; }
        (ChunkCoord Coord, int CellIndex)? Home   { get; }
        (ChunkCoord Coord, int CellIndex)? Job    { get; }

        // Lazy-allocate-or-fetch the per-action scratch bucket. Cast to
        // T via the new() constraint so the agent doesn't need to know
        // about every concrete subtype.
        T GetActionData<T>() where T : ActionData, new();

        // Display name for debug HUDs.
        string DisplayName { get; }
    }
}
