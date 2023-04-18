using UnityEngine;

public static class StaticResourcesLoader
{
    public static ComputeShader TriangleCandidateShader { get; private set; }

    [RuntimeInitializeOnLoadMethod(RuntimeInitializeLoadType.AfterAssembliesLoaded)]
    private static void LoadStaticAssets()
    {
        TriangleCandidateShader = Resources.Load<ComputeShader>("TriangleCandidate");
    }
}