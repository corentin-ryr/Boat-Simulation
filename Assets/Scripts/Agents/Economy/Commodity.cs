namespace Agents.Economy
{
    // Tradable goods carried by agents and stored in buildings. v1 keeps
    // this tiny — a richer set (Fish, Cloth, Wood, …) lives in the "later
    // commodities" phase.
    //
    // Money is intentionally NOT a Commodity: agents and buildings carry
    // money as a scalar `int Cash` field. Putting money in Inventory would
    // muddle "how much do I have to spend?" with "what am I carrying?",
    // and the price + payout math reads cleaner with cash separate.
    public enum Commodity : byte
    {
        Wheat = 0,
        Bread = 1,
    }

    // Sibling helper: enumeration + display names. Adding a Commodity =
    // one new enum entry + one new line in `All` + one new Name case.
    public static class Commodities
    {
        public static readonly Commodity[] All = { Commodity.Wheat, Commodity.Bread };
        public static int Count => All.Length;

        public static string Name(Commodity c) => c switch
        {
            Commodity.Wheat => "Wheat",
            Commodity.Bread => "Bread",
            _ => "?",
        };
    }
}
