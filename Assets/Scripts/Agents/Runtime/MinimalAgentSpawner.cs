using Agents.Core;
using UnityEngine;

namespace Agents.Runtime
{
    // Smallest viable spawner — drops `agentCount` agents at random
    // offsets from this transform with an EMPTY GoapBlueprint. The
    // agents register with AgentRegistry, decay their needs against
    // GameClock, and show up in the AgentDebugHud (F3). No goals, no
    // actions — pure runtime / registry / HUD validation for Phase E.
    //
    // Phase F replaces this with VillageSpawner, which instantiates the
    // real blueprint (sensors, actions, goals) plus assigns Home / Job
    // tiles to the spawned agents.
    public class MinimalAgentSpawner : MonoBehaviour
    {
        [Tooltip("How many smoke-test agents to spawn at Start.")]
        public int agentCount = 1;

        [Tooltip("Prefab carrying the Agent component (plus a visual mesh). " +
                 "If null, a fallback capsule is created.")]
        public GameObject agentPrefab;

        [Tooltip("Half-extents of the random XZ scatter around this transform.")]
        public Vector3 spawnArea = new Vector3(3f, 0f, 3f);

        void Start()
        {
            // Empty blueprint — runtime ticks Needs and registers with the
            // HUD, but no goal selection / planning runs.
            GoapBlueprint empty = new GoapBlueprint(null, null, null);

            for (int i = 0; i < agentCount; i++)
            {
                Vector3 pos = transform.position + new Vector3(
                    Random.Range(-spawnArea.x, spawnArea.x), 0f,
                    Random.Range(-spawnArea.z, spawnArea.z));

                GameObject go;
                if (agentPrefab != null)
                    go = Instantiate(agentPrefab, pos, Quaternion.identity);
                else
                    go = CreateFallbackVisual(pos, i);

                Agent agent = go.GetComponent<Agent>();
                if (agent == null) agent = go.AddComponent<Agent>();
                agent.Initialize(empty, $"Agent{i}", null, null);
            }
        }

        static GameObject CreateFallbackVisual(Vector3 pos, int idx)
        {
            GameObject go = GameObject.CreatePrimitive(PrimitiveType.Capsule);
            go.name = $"Agent{idx}";
            go.transform.position = pos;
            // Strip the capsule collider so the agent doesn't fight the
            // terrain mesh's raycast layer in TilePicker.
            Collider col = go.GetComponent<Collider>();
            if (col != null) Destroy(col);
            return go;
        }
    }
}
