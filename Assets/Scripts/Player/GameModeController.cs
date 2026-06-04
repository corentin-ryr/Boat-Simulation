using UnityEngine;
using TerrainGrid;

namespace Player
{
    // Single-source-of-truth for which control mode is active. Owns the Tab toggle
    // between gameplay (first/third-person + walkable avatar) and RTS (top-down
    // map editor) modes, and reconciles every subsystem that cares:
    //
    //   • Cameras: exactly one of {PlayerCamera GO, RtsCameraController GO} is
    //     active at any time. Activating/deactivating the parent GameObject
    //     disables the Camera component (so URP only renders the live one) AND
    //     the controller script in one move.
    //   • PlayerController: enabled in gameplay only — disabled in RTS so WASD
    //     doesn't double-drive the player.
    //   • TileSelector + TilePalette: enabled in RTS only. (TilePalette could
    //     stay enabled across modes for the persistent HUD, but disabling it in
    //     gameplay also clears its key-binding handlers and the OnGUI palette
    //     — cleaner.) TilePainter (gameplay-mode left-click paint) is enabled
    //     in gameplay only.
    //   • Cursor: locked + invisible in gameplay (existing PlayerController
    //     behaviour); free + visible in RTS so the user can click palette buttons
    //     and target tiles.
    //   • ChunkManager.cameraTarget: re-pointed to whichever mode's anchor drives
    //     streaming. The downstream Update() in ChunkManager already detects the
    //     chunk-coord change next frame and reissues the desired set.
    //   • Hover tint: pushes _HoverColor every mode flip, clears _HoverPrimalIdx
    //     on exit from RTS so a stale tint doesn't survive into gameplay. The
    //     per-cursor (_HoverPrimalIdx, _HoverChunkQR) updates happen on
    //     TileSelector each frame while RTS is active.
    //   • Outlines: flips ChunkSurface.ShowOutlines and calls
    //     RequestOutlineRefresh so every visible chunk's wedge-perimeter line
    //     mesh toggles on the same frame as Tab.
    //
    // Execution order: runs AFTER PlayerSpawner (which is at 100). The spawner
    // instantiates the player prefab at Start and writes the freshly-spawned
    // PlayerController + PlayerCamera GameObject + streaming anchor onto this
    // controller's Inspector fields; we then read them in our own Start to seed
    // initialMode. Running at 200 means PlayerSpawner is guaranteed to have
    // already filled those references in.
    [DefaultExecutionOrder(200)]
    public class GameModeController : MonoBehaviour
    {
        public enum Mode { Gameplay, Rts }

        [Header("Streaming")]
        public ChunkManager chunkManager;

        [Header("Gameplay")]
        [Tooltip("The avatar's PlayerController. Disabled in RTS so WASD doesn't double-up.")]
        public PlayerController playerController;
        [Tooltip("The gameplay camera GameObject. SetActive(false) in RTS so URP only " +
                 "renders the RTS camera.")]
        public GameObject playerCameraGo;
        [Tooltip("Anchor the streamer follows while in gameplay (typically the player's " +
                 "CameraAnchor). Set automatically by PlayerSpawner if left null at Start.")]
        public Transform gameplayStreamingAnchor;
        [Tooltip("The gameplay-mode TilePainter (left-click while walking). Disabled in RTS.")]
        public TilePainter tilePainter;

        [Header("RTS")]
        [Tooltip("The RTS camera's controller. Its GameObject is toggled active/inactive " +
                 "in lockstep with the gameplay camera.")]
        public RtsCameraController rtsCamera;
        [Tooltip("The RTS Pivot Transform — what streaming follows while in RTS. " +
                 "Defaults to rtsCamera.pivot if left null.")]
        public Transform rtsStreamingAnchor;
        [Tooltip("RTS-mode picker + hover. Disabled in gameplay.")]
        public TileSelector tileSelector;
        [Tooltip("Shared palette UI. Disabled in gameplay (no HUD, no key handlers).")]
        public TilePalette tilePalette;

        [Header("Toggle")]
        public KeyCode toggleKey = KeyCode.Tab;
        [Tooltip("Initial mode at scene start.")]
        public Mode initialMode = Mode.Gameplay;

        public Mode CurrentMode { get; private set; } = Mode.Gameplay;

        void Start()
        {
            // Fall-backs so the scene wires the minimum.
            if (rtsStreamingAnchor == null && rtsCamera != null) rtsStreamingAnchor = rtsCamera.pivot;
            if (gameplayStreamingAnchor == null && chunkManager != null) gameplayStreamingAnchor = chunkManager.cameraTarget;

            // Seed the RTS pivot at the gameplay anchor so a Tab press doesn't snap the
            // camera to wherever the RTS pivot happened to be parented in the prefab.
            // Done once at Start; subsequent toggles preserve the user's last position.
            if (rtsStreamingAnchor != null && gameplayStreamingAnchor != null)
                rtsStreamingAnchor.position = gameplayStreamingAnchor.position;

            // Apply initialMode unconditionally so all the gates land in a known state
            // (we can't trust whatever the prefab had toggled in the editor).
            ApplyMode(initialMode, instant: true);
        }

        void Update()
        {
            if (Input.GetKeyDown(toggleKey))
                ApplyMode(CurrentMode == Mode.Gameplay ? Mode.Rts : Mode.Gameplay, instant: false);
        }

        void ApplyMode(Mode mode, bool instant)
        {
            CurrentMode = mode;
            bool rts = mode == Mode.Rts;

            // Cameras: SetActive on the parent GameObject toggles the Camera
            // component, the controller script, and any children (like the
            // PlayerCamera's audio listener) in one move.
            if (playerCameraGo != null) playerCameraGo.SetActive(!rts);
            if (rtsCamera != null) rtsCamera.gameObject.SetActive(rts);

            // Pivot is auto-created in RtsCameraController.Awake. If the user
            // authored the RtsCamera GO as inactive, Awake hadn't run when our
            // Start fired — so rtsStreamingAnchor may still be null. Resolve it
            // now (we just activated the GO above, so Awake has run by this
            // point) and, when entering RTS, snap the pivot to the player's
            // current anchor. Two wins: bulletproof against edit-time inactive
            // RtsCamera, and Tab→RTS always centres on where the avatar is right
            // now (not where they spawned).
            if (rts && rtsCamera != null)
            {
                if (rtsStreamingAnchor == null) rtsStreamingAnchor = rtsCamera.pivot;
                if (rtsStreamingAnchor != null && gameplayStreamingAnchor != null)
                    rtsStreamingAnchor.position = gameplayStreamingAnchor.position;
            }

            // Locomotion: gameplay only. PlayerController lives on the spawned
            // player GameObject (which has to stay active for the avatar to exist
            // in the world even when not driven), so toggle the component, not the
            // GO.
            if (playerController != null) playerController.enabled = !rts;

            // Painters / palette: flip the whole GameObject's active state — that
            // covers both "user disabled the GO in the Inspector" and "user only
            // disabled the component". The component's .enabled stays at whatever
            // it was authored as; we don't fight it. SetActive is a no-op when
            // already in the requested state, so the per-frame Tab cost is nil.
            SetGoActive(tilePainter,   !rts);
            SetGoActive(tileSelector,  rts);
            SetGoActive(tilePalette,   rts);

            // Push the active camera onto each painter so they raycast through
            // the camera the user actually sees. Without this, the painter's
            // `cam` field gets cached at Awake (from Camera.main) and points at
            // whichever camera was main at scene start — usually the gameplay
            // one — which then gets SetActive(false) on the first Tab and
            // produces a stale frustum. Doing it from here in lockstep with the
            // mode flip means TileSelector / TilePainter don't need any
            // Inspector wiring for `cam` at all.
            Camera gameplayCam = playerCameraGo != null ? playerCameraGo.GetComponent<Camera>() : null;
            Camera rtsCam      = rtsCamera != null ? rtsCamera.GetComponent<Camera>() : null;
            if (tileSelector != null) tileSelector.cam = rtsCam;
            if (tilePainter  != null) tilePainter.cam  = gameplayCam;

            // Cursor: free in RTS, locked + hidden in gameplay (matches existing
            // PlayerController.Start behaviour, which we no longer let run on
            // gameplay re-entry since the controller is already enabled).
            Cursor.lockState = rts ? CursorLockMode.None : CursorLockMode.Locked;
            Cursor.visible   = rts;

            // Streaming anchor: re-point ChunkManager.cameraTarget. Its Update() will
            // re-evaluate the desired set on the next frame's WorldToChunk(cameraTarget).
            if (chunkManager != null)
            {
                Transform target = rts ? rtsStreamingAnchor : gameplayStreamingAnchor;
                if (target != null) chunkManager.cameraTarget = target;
            }

            // Hover lives in TerrainElevation.shader; push the colour + clear the
            // index on exit. _HoverPrimalIdx / _HoverChunkQR are driven by
            // TileSelector per cursor frame while RTS is active. Clearing the
            // index on exit from RTS avoids a stale tint surviving into gameplay
            // (where TileSelector is disabled and won't push -1 itself unless it
            // ran at least once).
            if (chunkManager != null)
                Shader.SetGlobalColor(HoverColorId, chunkManager.hoverColor);
            if (!rts) Shader.SetGlobalFloat(HoverPrimalIdxId, -1f);

            // Outline overlay lives on a per-chunk line-topology Mesh inside
            // ChunkSurface (BuildOutlineMesh). The flag is read every Apply for
            // newly-spawned presences; RequestOutlineRefresh pushes it into
            // already-live presences on the same frame as the Tab press.
            ChunkSurface.ShowOutlines = rts;
            chunkManager?.Surface?.RequestOutlineRefresh();

            _ = instant; // currently no animation; flag reserved for future tweens.
        }

        static readonly int HoverColorId     = Shader.PropertyToID("_HoverColor");
        static readonly int HoverPrimalIdxId = Shader.PropertyToID("_HoverPrimalIdx");

        // Tolerant of null Behaviours so the controller works even with optional
        // wiring left blank. Also tolerant of null GameObjects (a Behaviour without
        // a GO can't really happen, but the null-check keeps the call site clean).
        static void SetGoActive(Behaviour b, bool active)
        {
            if (b == null || b.gameObject == null) return;
            if (b.gameObject.activeSelf != active) b.gameObject.SetActive(active);
        }
    }
}
