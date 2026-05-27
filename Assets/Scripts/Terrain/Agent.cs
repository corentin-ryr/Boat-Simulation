using UnityEngine;

namespace TerrainGrid
{
    // Plain agent data — deliberately a struct in a flat array, not a MonoBehaviour per NPC, so
    // the store stays cache-friendly and scales to large crowds. Movement is continuous: the
    // agent steers toward Target and retargets on arrival. Rendering is a separate concern (only
    // near-camera agents are drawn), so an agent "exists" whether or not it's visible.
    public struct Agent
    {
        public Vector3 Position;
        public Vector3 Velocity;
        public Vector3 Target;   // current wander destination
        public float Speed;
    }
}
