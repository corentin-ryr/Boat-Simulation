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
        }
    }
}
