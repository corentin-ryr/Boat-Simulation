using UnityEngine;
using TerrainGrid;

namespace Player
{
    // Spawns the player at startup, places them above the sampled ground, and wires:
    //   • PlayerCamera.target/controller → the new instance
    //   • ChunkManager.cameraTarget       → the new instance's CameraAnchor
    //
    // Runs in Start with a higher DefaultExecutionOrder than ChunkManager so the streamer
    // has already primed the opening chunks and Elevation.Config is set (ChunkManager does
    // that in its Start). The user is expected to set ChunkManager.cameraTarget to this
    // spawner's Transform (or any Transform near the spawn point) in the inspector so the
    // *first* opening render set is built around the spawn point — we then hot-swap it to
    // the player after instantiation. Same set, no thrash.
    [DefaultExecutionOrder(100)]
    public class PlayerSpawner : MonoBehaviour
    {
        [Header("Refs")]
        public GameObject playerPrefab;
        public ChunkManager chunkManager;
        public PlayerCamera playerCamera;
        [Tooltip("Optional. If set, the spawner forwards the freshly-instantiated " +
                 "PlayerController and PlayerCamera's GameObject onto this mode " +
                 "controller — same pattern used to hand the camera anchor to " +
                 "ChunkManager. Leave null if the scene has no GameModeController.")]
        public GameModeController gameModeController;

        [Header("Spawn")]
        [Tooltip("Optional override. If null, spawns at this GameObject's position.")]
        public Transform spawnPoint;
        [Tooltip("Spawn this many metres above the sampled ground so the CharacterController " +
                 "settles cleanly on the first frame.")]
        public float spawnYBias = 0.5f;

        void Start()
        {
            if (playerPrefab == null)
            {
                Debug.LogWarning("PlayerSpawner: no playerPrefab assigned.", this);
                return;
            }

            Vector3 origin = spawnPoint != null ? spawnPoint.position : transform.position;
            float groundY = Elevation.Sample(origin.x, origin.z);
            Vector3 spawn = new Vector3(origin.x, groundY + spawnYBias, origin.z);

            GameObject p = Instantiate(playerPrefab, spawn, Quaternion.identity);
            PlayerController pc = p.GetComponent<PlayerController>();
            if (pc == null)
            {
                Debug.LogError("PlayerSpawner: playerPrefab is missing PlayerController.", this);
                return;
            }

            Transform anchor = pc.cameraAnchor != null ? pc.cameraAnchor : p.transform;

            if (playerCamera != null)
            {
                playerCamera.target = anchor;
                playerCamera.controller = pc;
            }
            if (chunkManager != null)
                chunkManager.cameraTarget = anchor;

            // Hand the freshly-spawned avatar + its camera to the mode controller
            // so Tab-toggling can disable them in RTS mode. Mirrors the existing
            // chunkManager.cameraTarget hand-off above: spawn-time wiring of a
            // scene-level controller that needs scene-instance references it can't
            // get from the prefab. The spawner also captures the gameplay
            // streaming anchor (the CameraAnchor we just gave ChunkManager) so the
            // controller can swap back to it on a Tab press out of RTS mode.
            if (gameModeController != null)
            {
                gameModeController.playerController = pc;
                if (playerCamera != null)
                    gameModeController.playerCameraGo = playerCamera.gameObject;
                gameModeController.gameplayStreamingAnchor = anchor;
            }
        }
    }
}
