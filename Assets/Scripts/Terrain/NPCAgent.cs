using UnityEngine;

namespace TerrainGrid
{
    // A physics-driven NPC. Lives on a prefab assigned to SimulationManager; one instance per
    // active NPC in the scene. Self-contained: holds its current wander target and steers its
    // own Rigidbody toward it. SimulationManager owns the crowd-level decisions (which targets
    // to pick, when to despawn) but never touches the Rigidbody directly.
    //
    // Steering runs in FixedUpdate so forces integrate cleanly with PhysX. Target updates from
    // the manager happen in Update — that's fine, the agent just reads the latest value next
    // physics step.
    [RequireComponent(typeof(Rigidbody))]
    public class NPCAgent : MonoBehaviour
    {
        [Tooltip("World-units per second the agent tries to reach.")]
        public float speed = 3f;

        [Tooltip("Cap on the steering force magnitude (Newtons). Prevents jittery overshoot.")]
        public float maxForce = 30f;

        // Wander destination in world space. Set by SimulationManager.
        public Vector3 Target;

        Rigidbody body;

        void Awake()
        {
            body = GetComponent<Rigidbody>();
        }

        void FixedUpdate()
        {
            float dt = Time.fixedDeltaTime;

            // Horizontal-only steering: gravity handles Y, the agent only decides where to walk.
            Vector3 here = transform.position;
            Vector3 to = Target - here; to.y = 0f;
            float dist = to.magnitude;
            if (dist < 1e-4f) return;

            Vector3 desiredVel = (to / dist) * speed;
            Vector3 currentVel = body.linearVelocity; currentVel.y = 0f;
            Vector3 deltaVel = desiredVel - currentVel;

            // Force needed to reach desiredVel in one step, clamped so a far target doesn't
            // launch the agent across the map.
            Vector3 force = body.mass * deltaVel / dt;
            float fmag = force.magnitude;
            if (fmag > maxForce) force *= maxForce / fmag;

            body.AddForce(force, ForceMode.Force);

            // Face the direction we're moving (visual only; physics doesn't care).
            Vector3 vel = body.linearVelocity; vel.y = 0f;
            if (vel.sqrMagnitude > 0.01f)
            {
                Quaternion look = Quaternion.LookRotation(vel.normalized, Vector3.up);
                body.MoveRotation(Quaternion.Slerp(body.rotation, look, 10f * dt));
            }
        }

        public bool ReachedTarget(float threshold)
        {
            Vector3 to = Target - transform.position; to.y = 0f;
            return to.sqrMagnitude < threshold * threshold;
        }
    }
}
