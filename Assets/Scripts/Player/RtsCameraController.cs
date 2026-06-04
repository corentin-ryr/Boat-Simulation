using UnityEngine;
using TerrainGrid;

namespace Player
{
    // Top-down RTS-style camera. Attached to a Camera GameObject that is the sibling
    // of the gameplay camera (PlayerCamera); enabled exclusively by GameModeController
    // when the player flips to RTS mode.
    //
    // Geometry: the camera orbits a logical Pivot point at fixed pitch (downward-
    // looking) and fixed yaw, panning horizontally and adjusting height via scroll
    // wheel. Pivot Y tracks the terrain via Elevation.Sample so the camera floats at
    // a consistent above-ground height when WASDing across hills.
    //
    // The Pivot Transform is what ChunkManager.cameraTarget points at while RTS is
    // active. WorldToChunk(pivot.position) drives streaming exactly the same way the
    // gameplay anchor did — no special-casing needed downstream.
    //
    // Controls (active only while this MonoBehaviour is enabled):
    //   WASD              pan pivot in camera-yaw-aligned XZ
    //   Scroll wheel      zoom (height) between minHeight and maxHeight
    //   Q / E             rotate yaw around pivot
    //   Shift             pan / zoom faster
    public class RtsCameraController : MonoBehaviour
    {
        [Header("Pivot")]
        [Tooltip("The look-at point the camera orbits. Created automatically as a child " +
                 "of this GameObject if left null.")]
        public Transform pivot;

        [Header("Pan")]
        [Tooltip("Pivot pan speed in m/s at minHeight (zoomed in). Linearly faster as the " +
                 "camera zooms out so high-altitude panning still feels responsive.")]
        public float panSpeed = 20f;
        [Tooltip("Multiplier on panSpeed at maxHeight relative to minHeight. 2 means twice " +
                 "as fast zoomed all the way out.")]
        public float panSpeedZoomFactor = 2f;
        [Tooltip("Shift-key multiplier on pan and zoom speeds.")]
        public float sprintMultiplier = 2.5f;

        [Header("Zoom")]
        public float minHeight = 8f;
        public float maxHeight = 80f;
        public float zoomSpeed = 8f;
        [Tooltip("Initial height of the camera above the pivot.")]
        public float initialHeight = 30f;

        [Header("Look")]
        [Tooltip("Pitch of the camera, degrees below horizontal. 90 = straight down, " +
                 "60 = a comfortable Townscaper-ish bird's-eye.")]
        public float pitchDegrees = 60f;
        [Tooltip("Yaw degrees per second when Q or E is held.")]
        public float yawDegreesPerSecond = 90f;
        [Tooltip("Initial yaw in degrees.")]
        public float initialYaw = 0f;

        [Header("Pivot ground tracking")]
        [Tooltip("Sample Elevation under the pivot every frame so the pivot Y tracks the " +
                 "terrain. Off → pivot Y stays at whatever it was initialised to (flat " +
                 "scrolling at fixed altitude, less smooth on hills).")]
        public bool trackTerrainHeight = true;

        float height;
        float yaw;

        void Awake()
        {
            if (pivot == null)
            {
                GameObject pivotGo = new GameObject("RtsPivot");
                pivotGo.transform.SetParent(transform.parent, false);
                pivotGo.transform.position = transform.position;
                pivot = pivotGo.transform;
            }
            height = Mathf.Clamp(initialHeight, minHeight, maxHeight);
            yaw = initialYaw;
            ApplyCameraTransform();
        }

        void OnEnable()
        {
            // When the mode controller flips us on, re-snap so the camera doesn't pop in
            // from a stale position (e.g. after a Time.timeScale = 0 / unfocused window).
            ApplyCameraTransform();
        }

        void Update()
        {
            float dt = Time.unscaledDeltaTime;
            float sprint = Input.GetKey(KeyCode.LeftShift) ? sprintMultiplier : 1f;

            // --- Zoom (scroll wheel) ---
            float wheel = Input.mouseScrollDelta.y;
            if (Mathf.Abs(wheel) > 0.001f)
            {
                height = Mathf.Clamp(height - wheel * zoomSpeed * sprint, minHeight, maxHeight);
            }

            // --- Yaw (Q / E) ---
            float yawInput = (Input.GetKey(KeyCode.E) ? 1f : 0f)
                           - (Input.GetKey(KeyCode.Q) ? 1f : 0f);
            yaw += yawInput * yawDegreesPerSecond * dt;

            // --- Pan (WASD), camera-yaw-aligned ---
            float h = Input.GetAxisRaw("Horizontal");
            float v = Input.GetAxisRaw("Vertical");
            if (Mathf.Abs(h) > 0.001f || Mathf.Abs(v) > 0.001f)
            {
                // Forward direction projected to XZ at the camera's current yaw.
                Quaternion yawRot = Quaternion.Euler(0f, yaw, 0f);
                Vector3 fwd = yawRot * Vector3.forward;
                Vector3 right = yawRot * Vector3.right;
                fwd.y = 0f; fwd.Normalize();
                right.y = 0f; right.Normalize();

                // Linear interpolation between panSpeed (at minHeight) and
                // panSpeed * panSpeedZoomFactor (at maxHeight).
                float t = Mathf.InverseLerp(minHeight, maxHeight, height);
                float zoomScale = Mathf.Lerp(1f, panSpeedZoomFactor, t);
                float speed = panSpeed * sprint * zoomScale;

                Vector3 wish = (fwd * v + right * h);
                if (wish.sqrMagnitude > 1f) wish.Normalize();

                Vector3 newPos = pivot.position + wish * speed * dt;
                pivot.position = newPos;
            }

            // Pivot tracks terrain Y so the camera holds a constant above-ground altitude
            // when panning over hills. Done after the pan so the height query uses the
            // freshly-moved XZ.
            if (trackTerrainHeight)
            {
                Vector3 pp = pivot.position;
                pp.y = Elevation.Sample(pp.x, pp.z);
                pivot.position = pp;
            }

            ApplyCameraTransform();
        }

        void ApplyCameraTransform()
        {
            // Place the camera at pivot + offset, where offset is built from yaw, pitch
            // and height. Computing it explicitly (rather than via LookAt on a parented
            // arm) keeps the math obvious and avoids accumulating-float drift in deep
            // hierarchies.
            Quaternion rot = Quaternion.Euler(pitchDegrees, yaw, 0f);
            // Direction from pivot to camera in the camera's local frame: backward and up.
            Vector3 offsetLocal = new Vector3(0f, 0f, -height / Mathf.Tan(pitchDegrees * Mathf.Deg2Rad));
            // The above puts the camera at horizontal distance = height/tan(pitch) along
            // -Z from the pivot in the camera's yaw plane, with `height` worth of vertical
            // separation. Equivalent to a polar offset (height, horizontalDist) rotated
            // by yaw and pointed at the pivot at the configured pitch.
            Vector3 horiz = Quaternion.Euler(0f, yaw, 0f) * offsetLocal;
            Vector3 camPos = pivot.position + horiz + Vector3.up * height;

            transform.SetPositionAndRotation(camPos, rot);
        }
    }
}
