using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public interface IWater {

    public abstract float GetWaterHeight(Vector3 position);

    Transform transform { get; }

}