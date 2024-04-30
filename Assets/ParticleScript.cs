using System.Collections;
using System.Collections.Generic;
using static System.Math;
using static System.Random;
using UnityEngine;
using Unity.Mathematics;
using static UnityEngine.Random;

public class ParticleScript : MonoBehaviour
{
    public int segments = 50;
    public float radius = 1f;

    LineRenderer line;

    public void Start()
    {
        line = GetComponent<LineRenderer>();
        if (line == null){
            line = gameObject.AddComponent<LineRenderer>();
            line.positionCount = segments + 1;
            line.useWorldSpace = false;
            line.startWidth = 0.1f;
            line.endWidth = 0.1f;
        }
        setRadius(0.2f);
        DrawCircle(new float2(0f, 0f));
    }

    public void DrawCircle(float2 curr_position)
    {
        float angle = 360f / segments;
        float theta = 0f;

        for (int i = 0; i < segments + 1; i++)
        {
            float x = Mathf.Sin(Mathf.Deg2Rad * theta) * radius;
            float y = Mathf.Cos(Mathf.Deg2Rad * theta) * radius;
            line.SetPosition(i, new Vector2(x, y) + new Vector2(curr_position.x, curr_position.y));
            theta += angle;
        }
        transform.position = new Vector2(curr_position.x, curr_position.y);
    }
    void setRadius(float new_radius){
         radius = new_radius;
    }
    
}
