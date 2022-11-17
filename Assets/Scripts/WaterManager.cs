using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;


public class WaterManager : MonoBehaviour
{

    List<Floater> floatingObjects = new List<Floater>();

    IWater water;

    // Start is called before the first frame update
    void Start()
    {
        water = FindObjectsOfType<MonoBehaviour>().OfType<IWater>().First();
    }

    
    public float GetWaterHeight(Vector3 position) {
        return water.GetWaterHeight(water.transform.InverseTransformDirection(position)) + water.transform.position.y;
    }


}
