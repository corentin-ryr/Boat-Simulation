using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Events;

public class EventManager 
{
    public static event UnityAction FinishedGeneratingMesh;
    public static void OnFinishedGeneratingMesh() => FinishedGeneratingMesh?.Invoke();
}

