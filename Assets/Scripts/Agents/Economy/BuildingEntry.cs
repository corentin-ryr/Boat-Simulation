using System.Collections.Generic;
using TerrainGrid;
using UnityEngine;

namespace Agents.Economy
{
    // Per-(chunk, cell) economic state. One BuildingEntry per role-bearing
    // tile in the world, owned and lifecycle-managed by BuildingRegistry.
    // Survives chunk eviction — the registry keeps the entry alive even
    // when the chunk it lives in unloads, so agents that were targeting
    // it don't lose their reference.
    //
    // What lives here:
    //   • Identity        — Coord, CellIndex, Role
    //   • Geometry        — WorldPos snapshot for cheap distance queries
    //   • Slot reservation — SlotsTotal + which agents have claimed them
    //   • Storage         — Stock (Field harvest, Bakery wheat/bread, etc.)
    //   • Trading         — Market (null for non-Market roles)
    //
    // What does NOT live here:
    //   • Mesh / GameObject — owned by ChunkSurface, lifecycle-detached.
    //   • TileProperty.Payload — gameplay metadata for the rendering layer;
    //     the BuildingEntry is the agent layer's source of truth.
    public class BuildingEntry
    {
        public readonly ChunkCoord Coord;
        public readonly int CellIndex;
        public BuildingRole Role { get; internal set; }
        public Vector3 WorldPos { get; internal set; }

        public int SlotsTotal { get; internal set; }

        // Claimed by agents at the start of an action that physically
        // occupies this tile (HarvestWheat, BakeBread, SellAtMarket).
        // Released on action Complete or Stop. Two-step ownership so a
        // dead agent doesn't permanently squat a slot.
        public readonly List<object> SlotsClaimed = new List<object>();

        // Resource storage. Bakery: stores Wheat (input) and Bread (output).
        // Field: holds harvested Wheat awaiting collection by a farmer.
        // Market: empty (use the MarketStockpile's own Stock instead).
        // House: empty in v1.
        public readonly Inventory Stock = new Inventory();

        // Market-only trading state. Null for every other role. Initialized
        // lazily by BuildingRegistry when a Market entry is created.
        public MarketStockpile Market;

        public int SlotsFree => SlotsTotal - SlotsClaimed.Count;
        public bool HasFreeSlot => SlotsFree > 0;

        public BuildingEntry(ChunkCoord coord, int cellIndex, BuildingRole role,
                             Vector3 worldPos, int slotsTotal)
        {
            Coord = coord;
            CellIndex = cellIndex;
            Role = role;
            WorldPos = worldPos;
            SlotsTotal = slotsTotal;
        }

        public override string ToString()
            => $"{BuildingRoles.Name(Role)}@{Coord}:{CellIndex} slots={SlotsClaimed.Count}/{SlotsTotal}";
    }
}
