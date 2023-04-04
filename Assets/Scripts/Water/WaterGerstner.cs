using System.Collections;
using System.Collections.Generic;
using UnityEngine;



[ExecuteInEditMode]
[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class WaterGerstner : MonoBehaviour, IWater
{
    MeshFilter meshFilter;
    Mesh sharedMeshCollider;
    public Mesh Mesh { get => meshFilter.mesh; }

    [Header("Parameters")]
    public int numGrid = 20;
    public float mapSize = 10;
    

    // Start is called before the first frame update
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();

        (Vector3[] vertices, int[] triangles, Vector2[] uvs) = MeshHelper.GenerateGridMesh(mapSize, numGrid, numGrid);
        meshFilter.mesh.vertices = vertices;
        meshFilter.mesh.triangles = triangles;
        meshFilter.mesh.uv = uvs;
        meshFilter.mesh.RecalculateNormals();

        gameObject.layer = 4; //4 is Water
    }


    private static Vector3 GerstnerWave(Vector3 position, float steepness, float wavelength, float speed, float direction)
    {
        direction = direction * 2 - 1;
        Vector2 d = new Vector2(Mathf.Cos(Mathf.PI * direction), Mathf.Sin(Mathf.PI * direction)).normalized;
        float k = 2 * Mathf.PI / wavelength;
        float a = steepness / k;
        float f = k * (Vector2.Dot(d, new Vector2(position.x, position.z)) - speed * Time.time);

        return new Vector3(d.x * a * Mathf.Cos(f), a * Mathf.Sin(f), d.y * a * Mathf.Cos(f));
    }

    public static Vector3 GetWaveDisplacement(Vector3 position, float steepness, float wavelength, float speed, float[] directions)
    {
        Vector3 offset = Vector3.zero;

        offset += GerstnerWave(position, steepness, wavelength, speed, directions[0]);
        offset += GerstnerWave(position, steepness, wavelength, speed, directions[1]);
        offset += GerstnerWave(position, steepness, wavelength, speed, directions[2]);
        offset += GerstnerWave(position, steepness, wavelength, speed, directions[3]);

        return offset;
    }

    public float GetWaterHeight(Vector3 position)
    {
        return GetWaveDisplacement(position, 0.1f, 1.5f, 0.5f, new float[] {0f, 1f, 0.6f, 0.3f}).y;
    }

}
