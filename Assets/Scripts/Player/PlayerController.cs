using UnityEngine;

namespace Player
{
    // Locomotion + body-yaw for the player avatar. Standard third-person split:
    //   • mouse X rotates THIS transform (the body), so WASD is body-forward-relative
    //   • mouse Y is owned by PlayerCamera (only the camera pitches; the body doesn't bend)
    //
    // The Mode enum gates whether on-foot input is applied. PlayerCamera flips it to FreeFly
    // (via F1) and back; while in FreeFly the controller goes idle and the camera does its
    // own movement.
    [RequireComponent(typeof(CharacterController))]
    public class PlayerController : MonoBehaviour
    {
        public enum Mode { OnFoot, FreeFly }

        [Header("Locomotion")]
        public float walkSpeed = 4f;
        public float sprintMultiplier = 1.6f;
        public float jumpHeight = 1.2f;
        public float gravity = -9.81f;
        [Tooltip("Small constant downward velocity while grounded — keeps the controller " +
                 "pressed onto slopes so isGrounded stays true on downhill strides.")]
        public float groundedSnapDown = 2f;

        [Header("Look")]
        public float mouseSensitivityX = 2f;

        [Header("Refs")]
        [Tooltip("Empty child at head height. The camera pivots around this point and the " +
                 "ChunkManager uses it as its streaming anchor.")]
        public Transform cameraAnchor;

        // Read by PlayerCamera; PlayerCamera writes when toggling free-fly.
        [HideInInspector] public Mode mode = Mode.OnFoot;

        CharacterController cc;
        Vector3 velocity;
        float yaw;

        void Awake()
        {
            cc = GetComponent<CharacterController>();
            yaw = transform.eulerAngles.y;
        }

        void Start()
        {
            // Lock + hide the cursor for mouselook. Press Esc in editor to unlock.
            Cursor.lockState = CursorLockMode.Locked;
            Cursor.visible = false;
        }

        void Update()
        {
            if (mode != Mode.OnFoot) return;

            // --- Yaw: mouse X rotates the body so movement is camera-aligned ---
            yaw += Input.GetAxisRaw("Mouse X") * mouseSensitivityX;
            transform.rotation = Quaternion.Euler(0f, yaw, 0f);

            // --- Horizontal movement (body-relative WASD) ---
            float h = Input.GetAxisRaw("Horizontal");
            float v = Input.GetAxisRaw("Vertical");
            Vector3 wish = transform.forward * v + transform.right * h;
            if (wish.sqrMagnitude > 1f) wish.Normalize();
            float speed = walkSpeed * (Input.GetKey(KeyCode.LeftShift) ? sprintMultiplier : 1f);

            // --- Gravity + jump ---
            if (cc.isGrounded)
            {
                velocity.y = -groundedSnapDown;
                if (Input.GetButtonDown("Jump"))
                    velocity.y = Mathf.Sqrt(jumpHeight * -2f * gravity);
            }
            else
            {
                velocity.y += gravity * Time.deltaTime;
            }

            Vector3 move = wish * speed + Vector3.up * velocity.y;
            cc.Move(move * Time.deltaTime);
        }
    }
}
