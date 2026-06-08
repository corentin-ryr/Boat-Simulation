using TerrainGrid;

namespace Agents.Economy
{
    // Economic role of a tile in the agent simulation. Derived from
    // TileKind: most kinds map 1-to-1 (Field → Field, Market → Market),
    // role-bearing sub-kinds (House, Bakery) inherit their own role,
    // non-economic kinds (Default, Road, Building generic, Lighthouse)
    // map to None.
    //
    // Adding an economic TileKind = one new enum entry here + the
    // TileKind → BuildingRole switch + the default-slots table.
    public enum BuildingRole : byte
    {
        None   = 0,
        House  = 1,  // sleep target, family residence
        Field  = 2,  // wheat production site
        Bakery = 3,  // wheat → bread workshop
        Market = 4,  // trading hub with stockpile + prices
        Dock   = 5,  // future: fish production / caravan terminus
    }

    public static class BuildingRoles
    {
        // Source-of-truth mapping. The reverse direction (role → TileKind)
        // isn't always meaningful — some roles could map to several sub-
        // kinds in the future (e.g. Workshop kind for the "Bakery" role).
        public static BuildingRole FromTileKind(TileKind kind) => kind switch
        {
            TileKind.House  => BuildingRole.House,
            TileKind.Field  => BuildingRole.Field,
            TileKind.Bakery => BuildingRole.Bakery,
            TileKind.Market => BuildingRole.Market,
            TileKind.Dock   => BuildingRole.Dock,
            _               => BuildingRole.None,
        };

        // How many agents can simultaneously occupy a working slot at this
        // role's tile. Drives the slot-claim system in BuildingRegistry —
        // two farmers can't claim the same Field, but two could claim
        // adjacent Fields. House slots = residents; Market slots =
        // simultaneous customers + traders standing at the counter.
        public static int DefaultSlots(BuildingRole role) => role switch
        {
            BuildingRole.House  => 1,
            BuildingRole.Field  => 2,
            BuildingRole.Bakery => 1,
            BuildingRole.Market => 4,
            BuildingRole.Dock   => 2,
            _                   => 0,
        };

        public static string Name(BuildingRole r) => r switch
        {
            BuildingRole.House  => "House",
            BuildingRole.Field  => "Field",
            BuildingRole.Bakery => "Bakery",
            BuildingRole.Market => "Market",
            BuildingRole.Dock   => "Dock",
            _                   => "None",
        };
    }
}
