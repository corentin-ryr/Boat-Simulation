using System;
using UnityEngine;
using World;

namespace Agents.Runtime
{
    // Hunger / Energy state for one agent, decayed against in-game time.
    // Wall-time decay is the wrong driver here — fast-forwarding the
    // GameClock (F2 cycle) should drive needs forward too, not just the
    // agents' visible motion. Tracks last seen TimeOfDay in game-hours
    // so a 16× timeScale produces 16× faster decay, no edge-case math.
    //
    // Decay rates are hourly (game-time hours). Hunger and Energy both
    // start at 100; Eat / Sleep bring them back up.
    [Serializable]
    public class AgentNeeds
    {
        [Range(0f, 100f)] public float Hunger = 100f;
        [Range(0f, 100f)] public float Energy = 100f;

        [Tooltip("Hunger lost per in-game hour while awake.")]
        public float HungerDecayPerHour = 8f;
        [Tooltip("Energy lost per in-game hour while awake.")]
        public float EnergyDecayPerHour = 4f;
        [Tooltip("Energy gained per in-game hour while sleeping.")]
        public float EnergyRecoveryPerHour = 30f;
        [Tooltip("Hunger gained per Bread eaten.")]
        public float HungerPerBread = 30f;

        // Tracks the last sampled total elapsed game-hours so cross-day
        // wrap doesn't trip up the diff. Negative means "uninitialized".
        float lastTotalHours = -1f;

        public void Tick(bool isSleeping)
        {
            GameClock clock = GameClock.Instance;
            if (clock == null) return;

            float nowTotal = clock.DayOfYear * GameClock.HoursPerDay + clock.TimeOfDayHours;
            if (lastTotalHours < 0f) { lastTotalHours = nowTotal; return; }

            float dHours = nowTotal - lastTotalHours;
            if (dHours <= 0f) return;
            lastTotalHours = nowTotal;

            // Hunger decays whether asleep or awake — you wake up hungry.
            Hunger = Mathf.Clamp(Hunger - HungerDecayPerHour * dHours, 0f, 100f);

            Energy = isSleeping
                ? Mathf.Clamp(Energy + EnergyRecoveryPerHour * dHours, 0f, 100f)
                : Mathf.Clamp(Energy - EnergyDecayPerHour * dHours, 0f, 100f);
        }

        public void EatBread() => Hunger = Mathf.Clamp(Hunger + HungerPerBread, 0f, 100f);

        public bool IsHungry(float threshold = 40f) => Hunger < threshold;
        public bool IsTired(float threshold = 35f) => Energy < threshold;
        public bool IsRested(float threshold = 90f) => Energy >= threshold;
        public bool IsFed(float threshold = 70f) => Hunger >= threshold;
    }
}
