using UnityEngine;

namespace Player
{
    // Third-person spring-arm orbit around the player's CameraAnchor, with a free-fly mode
    // toggled by F1.
    //
    // Yaw split: in third-person, this camera does NOT own yaw — the player body already
    // rotated to mouse X (in PlayerController), and we just read its eulerAngles.y. We own
    // pitch (mouse Y) and the arm length / occlusion test. In free-fly, this camera owns
    // both yaw and pitch independently of the player.
    public class PlayerCamera : MonoBehaviour
    {
        [Header("Refs")]
        [Tooltip("The player's CameraAnchor (head-height empty). Assigned by PlayerSpawner.")]
        public Transform target;
        [Tooltip("PlayerController whose mode this camera reads + flips on F1. " +
                 "Assigned by PlayerSpawner.")]
        public PlayerController controller;

        [Header("Third-person")]
        public float armLength = 4f;
        [Tooltip("Sphere-cast radius. Wider means the camera pulls in earlier near obstacles.")]
        public float armRadius = 0.3f;
        public float minPitch = -30f;
        public float maxPitch = 60f;
        public float mouseSensitivityY = 2f;
        [Tooltip("Which layers count as occluders. Default ~0 = everything. Set to your " +
                 "terrain layer if you have one.")]
        public LayerMask occlusionMask = ~0;

        [Header("Free-fly")]
        public float freeFlySpeed = 12f;
        public float freeFlyBoostMultiplier = 4f;
        public float freeFlyMouseSensitivity = 2f;

        float pitch;
        // Free-fly state — owns its own yaw + pitch so the view doesn't fight the player body.
        float freeYaw, freePitch;

        void LateUpdate()
        {
            if (Input.GetKeyDown(KeyCode.F1)) ToggleFreeFly();

            if (controller != null && controller.mode == PlayerController.Mode.FreeFly)
                UpdateFreeFly();
            else
                UpdateThirdPerson();
        }

        void UpdateThirdPerson()
        {
            if (target == null) return;

            pitch -= Input.GetAxisRaw("Mouse Y") * mouseSensitivityY;
            pitch = Mathf.Clamp(pitch, minPitch, maxPitch);

            // Yaw mirrors the player body (the player already aligned to mouse X this frame).
            Quaternion rot = Quaternion.Euler(pitch, target.eulerAngles.y, 0f);
            Vector3 desired = target.position - rot * Vector3.forward * armLength;

            // Sphere-cast from the target outward; if we hit something, pull the camera in.
            Vector3 dir = (desired - target.position).normalized;
            if (Physics.SphereCast(target.position, armRadius, dir, out RaycastHit hit,
                                   armLength, occlusionMask, QueryTriggerInteraction.Ignore))
                desired = target.position + dir * Mathf.Max(0f, hit.distance - 0.05f);

            transform.SetPositionAndRotation(desired, rot);
        }

        void UpdateFreeFly()
        {
            freeYaw   += Input.GetAxisRaw("Mouse X") * freeFlyMouseSensitivity;
            freePitch -= Input.GetAxisRaw("Mouse Y") * freeFlyMouseSensitivity;
            freePitch = Mathf.Clamp(freePitch, -89f, 89f);
            Quaternion rot = Quaternion.Euler(freePitch, freeYaw, 0f);

            float speed = freeFlySpeed * (Input.GetKey(KeyCode.LeftShift) ? freeFlyBoostMultiplier : 1f);
            Vector3 fwd = rot * Vector3.forward;
            Vector3 right = rot * Vector3.right;
            float h = Input.GetAxisRaw("Horizontal");
            float v = Input.GetAxisRaw("Vertical");
            float y = (Input.GetKey(KeyCode.E) ? 1f : 0f) - (Input.GetKey(KeyCode.Q) ? 1f : 0f);

            Vector3 delta = (fwd * v + right * h + Vector3.up * y) * speed * Time.deltaTime;
            transform.SetPositionAndRotation(transform.position + delta, rot);
        }

        void ToggleFreeFly()
        {
            if (controller == null) return;
            if (controller.mode == PlayerController.Mode.OnFoot)
            {
                controller.mode = PlayerController.Mode.FreeFly;
                // Seed free-fly orientation from current camera so the view doesn't snap.
                Vector3 e = transform.eulerAngles;
                freePitch = NormalizeAngle(e.x);
                freeYaw = e.y;
            }
            else
            {
                controller.mode = PlayerController.Mode.OnFoot;
            }
        }

        static float NormalizeAngle(float a) => a > 180f ? a - 360f : a;
    }
}
