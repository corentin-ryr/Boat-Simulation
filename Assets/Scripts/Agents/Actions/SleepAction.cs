using System.Collections.Generic;
using Agents.Core;
using Agents.Runtime;

namespace Agents.Actions
{
    // Lie down at Home and refill Energy. While sleeping, AgentNeeds.Tick
    // adds EnergyRecoveryPerHour × elapsed-game-hours to Energy instead
    // of decaying it. Action.Start flips Agent.IsSleeping = true so the
    // Tick reads "sleeping"; Action.Complete / Stop flips it back.
    //
    // Symbolic contract:
    //   Pre  : AtHome = true, Rested = false
    //   Eff  : Rested = true
    //
    // Real-world contract: Perform polls Energy and returns true once
    // it crosses the IsRested threshold (~90). Wall-time elapsed depends
    // on timeScale: at default 30 wall-s / game-hour and 30 energy / hour
    // recovery, going from Energy=0 to Energy=90 takes ~90 wall-seconds.
    public sealed class SleepAction : GoapAction
    {
        public override string Name => "Sleep";

        static readonly KeyValuePair<string, object>[] PreCond =
        {
            new KeyValuePair<string, object>(WorldKey.AtHome, true),
            new KeyValuePair<string, object>(WorldKey.Rested, false),
        };
        static readonly KeyValuePair<string, object>[] Eff =
        {
            new KeyValuePair<string, object>(WorldKey.Rested, true),
        };
        public override IEnumerable<KeyValuePair<string, object>> Preconditions => PreCond;
        public override IEnumerable<KeyValuePair<string, object>> Effects => Eff;
        public override float Cost(IGoapAgent agent) => 1f;

        public override void Start(IGoapAgent agent)
        {
            Agent concrete = agent as Agent;
            if (concrete != null) concrete.IsSleeping = true;
        }

        public override bool Perform(IGoapAgent agent, float dt)
        {
            Agent concrete = agent as Agent;
            if (concrete == null) return true;
            return concrete.Needs.IsRested();
        }

        public override void Complete(IGoapAgent agent)
        {
            Agent concrete = agent as Agent;
            if (concrete != null) concrete.IsSleeping = false;
        }

        public override void Stop(IGoapAgent agent)
        {
            Agent concrete = agent as Agent;
            if (concrete != null) concrete.IsSleeping = false;
        }
    }
}
