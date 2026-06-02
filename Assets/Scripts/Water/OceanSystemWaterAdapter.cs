using OceanSystem;
using UnityEngine;

// Bridges the gasgiant/Ocean-URP package's FFT ocean to the project's existing
// IWater interface so WaterManager / Floater / PressureFloater can query height
// against the new ocean without any changes to the buoyancy code.
//
// WaterManager.Start does FindObjectsOfType<MonoBehaviour>().OfType<IWater>().First();
// so just having this component on any active scene GameObject (we put it on the
// Ocean prefab) is enough for the wiring. The intended layout is to sit on the
// same GameObject as the OceanSimulation — RequireComponent enforces that and
// lets us auto-resolve the reference in Awake when it isn't set explicitly.
[RequireComponent(typeof(OceanSimulation))]
[AddComponentMenu("Ocean/Ocean System Water Adapter")]
public class OceanSystemWaterAdapter : MonoBehaviour, IWater
{
    [Tooltip("Source simulation. Auto-resolved from the same GameObject if left null.")]
    public OceanSimulation simulation;

    void Awake()
    {
        if (simulation == null) simulation = GetComponent<OceanSimulation>();
    }

    public float GetWaterHeight(Vector3 position)
    {
        // OceanSimulation.Collision is null until the first simulation step has
        // run, and the underlying readback textures need a couple of frames before
        // they hold valid data — OceanCollision.GetWaterHeight returns 0 internally
        // until then, but the Collision instance itself can also be null on the
        // very first frames after scene load. Guard both.
        if (simulation == null || simulation.Collision == null) return 0f;
        return simulation.Collision.GetWaterHeight(position);
    }
}
