using System;
using UnityEngine;

namespace World
{
    // Singleton game-time driver. The sole source of "what time is it?" for
    // gameplay code — agent needs decay, daily market price ticks, day-night
    // visibility, etc. all subscribe here instead of reading Time.deltaTime
    // directly so the simulation can be paused / fast-forwarded / replayed
    // from one switch.
    //
    // Model: one in-game day = HoursPerDay (24) in-game hours; one in-game
    // hour = secondsPerGameHour wall-clock seconds. Hour 0 starts at
    // midnight; the daytime window is [DaytimeStartHour, DaytimeEndHour).
    //
    // Events:
    //   OnHour(int hour)     — fired the moment GameHour transitions to the
    //                          new value (e.g. 7 → 8). Hour is in [0, 24).
    //                          Used for: needs decay, daily price tick (at
    //                          hour 0), shift changes.
    //   OnNewDay(int day)    — fired right after OnHour when GameHour wraps
    //                          past 24. day = DayOfYear after the wrap.
    //
    // Pause / fast-forward:
    //   timeScale = 0  → clock frozen.
    //   timeScale > 1  → wall-clock seconds compress; one F2 press cycles
    //                    1× → 4× → 16× → 1× for debug.
    public class GameClock : MonoBehaviour
    {
        public const int HoursPerDay = 24;
        public const int DaytimeStartHour = 6;   // 06:00 inclusive
        public const int DaytimeEndHour = 21;    // 21:00 exclusive

        public static GameClock Instance { get; private set; }

        [Header("Game time")]
        [Tooltip("Wall-clock seconds per in-game hour. 30 = a 24h day in 12 real minutes; " +
                 "60 = a 24h day in 24 real minutes.")]
        public float secondsPerGameHour = 30f;

        [Tooltip("Initial in-game hour at scene start. Default 7 = soon after dawn so the " +
                 "village wakes up immediately.")]
        [Range(0, HoursPerDay - 1)]
        public int startHour = 7;

        [Header("Debug")]
        [Tooltip("Multiplier applied to the clock advance. Press cycleSpeedKey to cycle 1× / 4× / 16×.")]
        public float timeScale = 1f;
        public KeyCode cycleSpeedKey = KeyCode.F2;

        [Tooltip("Show an on-screen HUD with the current time and day.")]
        public bool drawHud = true;

        // Fractional accumulator. Whole-hour transitions fire OnHour; whole-day
        // wraps fire OnNewDay. Subscribers should treat the integer hour /
        // day as the source of truth, not the float.
        float hourAccumulator;
        int gameHour;
        int dayOfYear;

        // Public read-only view of the in-game state.
        public int GameHour => gameHour;
        public int DayOfYear => dayOfYear;
        public float HourFraction => hourAccumulator;          // [0, 1)
        public bool IsDaytime => gameHour >= DaytimeStartHour && gameHour < DaytimeEndHour;
        public bool IsNighttime => !IsDaytime;
        // Continuous time-of-day for shader feeds / smooth interpolation.
        public float TimeOfDayHours => gameHour + hourAccumulator;

        public event Action<int> OnHour;     // hour after transition (0..23)
        public event Action<int> OnNewDay;   // day after wrap

        void Awake()
        {
            // Single-instance pattern. The clock is scene-bound; if a second
            // scene Additive-loads its own GameClock, the new one loses.
            if (Instance != null && Instance != this)
            {
                Destroy(this);
                return;
            }
            Instance = this;
            gameHour = Mathf.Clamp(startHour, 0, HoursPerDay - 1);
            hourAccumulator = 0f;
            dayOfYear = 0;
        }

        void OnDestroy()
        {
            if (Instance == this) Instance = null;
        }

        void Update()
        {
            if (Input.GetKeyDown(cycleSpeedKey)) CycleTimeScale();

            float secsPerHour = Mathf.Max(0.01f, secondsPerGameHour);
            float dHours = Time.deltaTime * Mathf.Max(0f, timeScale) / secsPerHour;
            if (dHours <= 0f) return;

            hourAccumulator += dHours;
            while (hourAccumulator >= 1f)
            {
                hourAccumulator -= 1f;
                gameHour++;
                if (gameHour >= HoursPerDay)
                {
                    gameHour = 0;
                    dayOfYear++;
                    OnHour?.Invoke(gameHour);
                    OnNewDay?.Invoke(dayOfYear);
                }
                else
                {
                    OnHour?.Invoke(gameHour);
                }
            }
        }

        void CycleTimeScale()
        {
            if (timeScale < 2f) timeScale = 4f;
            else if (timeScale < 8f) timeScale = 16f;
            else timeScale = 1f;
        }

        void OnGUI()
        {
            if (!drawHud) return;
            const float w = 220f;
            const float h = 48f;
            float x = Screen.width - w - 10f;
            float y = 10f;
            GUI.Box(new Rect(x, y, w, h), "");
            string period = IsDaytime ? "day" : "night";
            GUI.Label(new Rect(x + 8f, y + 4f, w - 16f, 20f),
                      $"Day {dayOfYear}  {gameHour:D2}:{(int)(hourAccumulator * 60f):D2}  ({period})");
            GUI.Label(new Rect(x + 8f, y + 24f, w - 16f, 20f),
                      $"speed {timeScale:0.#}×   [{cycleSpeedKey}]");
        }
    }
}
