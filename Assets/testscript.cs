using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using static System.Math;
using static System.Random;
using UnityEngine;
using Unity.Mathematics;
using static UnityEngine.Random;

public class testscript : MonoBehaviour
{
    public float gravity;
    float smoothingRadius;
    public static int num_particles = 500;
    ParticleScript[] particles;

    public float2[] positions;
    public float2[] velocities;

    public ArrayList[] gridArr;
    public int[] sector_list;
    public static int numSectors = 25;

    float[] densities;

    ParticleScript particle;
    public float particle_size;
    public float particle_mass = 1.5f;
    public float damping_factor;

    public float ideal_density;
    public float pressure_factor;

    float2 bounds_size = new float2(15f, 15f);

    // Start is called before the first frame update
    void Start()
    {
        gravity = 4f;
        ideal_density = 1.1f * num_particles / (bounds_size.x * bounds_size.y);
        pressure_factor = (5.3f * gravity) / ideal_density;
        smoothingRadius = 1f;
        particle_size = 1f;
        damping_factor = 0.35f;

        densities = new float[num_particles];
        positions = new float2[num_particles];
        velocities = new float2[num_particles];
        gridArr = new ArrayList[25];
        for (int grid_number = 0;grid_number<numSectors;grid_number++){
            gridArr[grid_number] = new ArrayList();
        }
        sector_list = new int[num_particles];

        particle = FindObjectOfType<ParticleScript>();
        if (particle!=null){
            particle.DrawCircle(new float2(0f, 0f));
        }
        else{
            Debug.Log("bad stuff");
        }
        particles = new ParticleScript[num_particles];
        particles[0] = particle;
        positions[0] = float2.zero;
        velocities[0] = float2.zero;


        for(int i = 1;i<num_particles;i++){

            float random_x = Range(-7f, 7f);
            float random_y = Range(-7f, 7f);

            float2 particle_pos = new float2(random_x, random_y);


            Object ith_particle = (ParticleScript) Instantiate(particle, new Vector3(random_x, random_y, 0f), Quaternion.identity);
            particles[i] = (ParticleScript) ith_particle;
            particles[i].Start();
            particles[i].DrawCircle(particle_pos);
            positions[i] = particle_pos;
            velocities[i] = float2.zero;

            int particle_sector = getGridSector(random_x, random_y);
            gridArr[particle_sector].Add(i);
            sector_list[i] = particle_sector;
        }
        for (int i = 0;i<num_particles;i++){
            densities[i] = getDensityNew(positions[i]);
        }
    }

    // Update is called once per frame
    void Update()
    {
        float time_step = Time.deltaTime;
        updateGridArray();

        Parallel.For(0,num_particles,i=>{
            velocities[i]+= (new float2(0f,-1f)) * gravity * time_step;
            densities[i] = getDensityNew(positions[i]);
        });

        for (int i = 0;i<num_particles;i++){
            float2 pressureForce = calculatePressureForceNew(i);
            float2 a_pressure = pressureForce /densities[i];
            velocities[i] += (a_pressure) * time_step;
        }
        /*
        for (int i = 0;i<num_particles;i++){
            velocities[i]+= (new float2(0f,-1f)) * gravity * time_step;
            densities[i] = getDensityNew(positions[i]);
        }
        for (int i = 0;i<num_particles;i++){
            positions[i]+= velocities[i] * time_step;
            FixCollisions(i);
        }
        */
        Parallel.For(0,num_particles,i=>{
            positions[i]+= velocities[i] * time_step;
            FixCollisions(i);
        });

        for (int i = 0;i<num_particles;i++){
            particles[i].DrawCircle(positions[i]);
        }

        

    }

    void FixCollisions(int i){
        float2 half_bound_size = bounds_size / 2 - (new float2(1,1)) * particle_size;

        if (Abs(positions[i].x) > half_bound_size.x){
            positions[i].x = half_bound_size.x * Sign(positions[i].x);
            velocities[i].x *=-1 * damping_factor;
        }
        if (Abs(positions[i].y) > half_bound_size.y){
            positions[i].y = half_bound_size.y * Sign(positions[i].y);
            velocities[i].y *=-1.5f * damping_factor;
        }
    }
    float smoothingKernelNew(float radius, float distance){

        float h = radius / 2f;

        float scaling_factor = (float) (1f / (PI * h * h * h));

        if (distance >=0 && distance <=h){
            return (float) (scaling_factor * (1 - 1.5*((distance / h) * (distance / h)) + 0.75*((distance/h) * (distance/h) * (distance/h))));
        }
        else if (distance > h && distance <= radius){
            return (float) (scaling_factor * (0.25*((2-(radius/h))* (2-(radius/h)) * (2-(radius/h)))));
        }
        else{
            return 0f;
        }
    }

    float smoothingKernelDerivativeNew(float radius, float distance){
        float h = radius / 2f;

        if (distance >=0 && distance <=h){
            return (float) (((9f / (4f * PI * Pow(h, 6))) * (distance * distance)) - ((6 * distance) /(2 * PI * Pow(h, 5))));
        }
        else if (distance > h && distance <= radius){
            float scaling_factor = (float) (1f / (4* PI * h* h * h));

            return scaling_factor * (((-3 * distance * distance)/(h * h * h)) + (12 * distance / (h * h)) - (12 / h));
        }
        return 0f;
    }
    float getDensityNew(float2 density_position){
        float x_position = density_position.x;
        float y_position = density_position.y;

        int current_sector = getGridSector(x_position, y_position);
        int[] nearby_sectors = getNearbySectors(current_sector);
        float density = 0f;
        foreach (int sector in nearby_sectors){
            ArrayList given_sector = gridArr[sector];

            foreach (int particle_index in given_sector){
                float distance = magnitude(positions[particle_index] - density_position);
                float influence = smoothingKernelNew(smoothingRadius, distance);

                density+= particle_mass * influence;
            }
        }

        return (float) Max(density, 0.01f);
    }

    float magnitude(float2 vector){
        return (float) Pow(Pow(vector.x, 2) + Pow(vector.y, 2), 0.5);
    }


    //3 x 3 sectors of area 9 each (25 such sectors)

    // 0  1  2  3  4
    // 5  6  7  8  9
    //10 11 12 13 14
    //15 16 17 18 19
    //20 21 22 23 24

    int getGridSector(float x_coord, float y_coord){
        int[,] sector_map = new int[5,5]{
            {0,1,2,3,4},
            {5,6,7,8,9},
            {10,11,12,13,14},
            {15,16,17,18,19},
            {20,21,22,23,24}
        };

        float x_val = x_coord+7.50f;
        float y_val = y_coord+7.50f;

        x_val = Min(x_val, 14.99f);
        y_val = Min(y_val, 14.99f);

        x_val/=3;
        y_val/=3;

        int x_sector = (int) x_val;
        int y_sector = (int) y_val;

        y_sector = 4-y_sector;

        return sector_map[y_sector,x_sector];

    }

    int[] getNearbySectors(int sector){
        int[] to_return;

        switch(sector){
            case 0:
                to_return = new int[] {0,1,5,6};
                break;
            case 1:
                to_return = new int[] {0,1,2,5,6,7};
                break;
            case 2:
                to_return = new int[] {1,2,3,6,7,8};
                break;
            case 3:
                to_return = new int[] {2,3,4,7,8,9};
                break;
            case 4:
                to_return = new int[] {3,4,8,9};
                break;
            case 5:
                to_return = new int[] {0,1,5,6,10,11};
                break;
            case 6:
                to_return = new int[] {0,1,2,5,6,7,10,11,12};
                break;
            case 7:
                to_return = new int[] {1,2,3,6,7,8,11,12,13};
                break;
            case 8:
                to_return = new int[] {2,3,4,7,8,9,12,13,14};
                break;
            case 9:
                to_return = new int[] {3,4,8,9,13,14};
                break;
            case 10:
                to_return = new int[] {5,6,10,11,15,16};
                break;
            case 11:
                to_return = new int[] {5,6,7,10,11,12,15,16,17};
                break;
            case 12:
                to_return = new int[] {6,7,8,11,12,13,16,17,18};
                break;
            case 13:
                to_return = new int[] {7,8,9,12,13,14,17,18,19};
                break;
            case 14:
                to_return = new int[] {8,9,13,14,18,19};
                break;
            case 15:
                to_return = new int[] {10,11,15,16,20,21};
                break;
            case 16:
                to_return = new int[] {10,11,12,15,16,17,20,21,22};
                break;
            case 17:
                to_return = new int[] {11,12,13,16,17,18,21,22,23};
                break;
            case 18:
                to_return = new int[] {12,13,14,17,18,19,22,23,24};
                break;
            case 19:
                to_return = new int[] {13,14,18,19,23,24};
                break;
            case 20:
                to_return = new int[] {15,16,20,21};
                break;
            case 21:
                to_return = new int[] {15,16,17,20,21,22};
                break;
            case 22:
                to_return = new int[] {16,17,18,21,22,23};
                break;
            case 23:
                to_return = new int[] {17,18,19,22,23,24};
                break; 
            case 24:
                to_return = new int[] {18,19,23,24};
                break;
            default:
                Debug.Log("Problems have occurred");
                to_return = new int[1];
                to_return[0] = 0;
                break;
        }

        return to_return;
    }

    void updateGridArray(){
        for (int i = 0;i<num_particles;i++){
            float x_position = positions[i].x;
            float y_position = positions[i].y;

            int new_sector = getGridSector(x_position, y_position);
            if (new_sector == sector_list[i]){
                continue;
            }
            
            gridArr[sector_list[i]].Remove(i);
            gridArr[new_sector].Add(i);
            sector_list[i] = new_sector;
        }
    }

    float arbitraryFunc(float2 pos){
        return Mathf.Cos(pos.y - 3 + Mathf.Sin(pos.x));
    }

    float2 calculatePressureForceNew(int particle_index){
        float2 pressureForce = new float2(0f, 0f);

        float2 particle_position = positions[particle_index];
        int particle_sector = getGridSector(particle_position.x, particle_position.y);
        int[] nearby_sectors = getNearbySectors(particle_sector);

        foreach (int sector in nearby_sectors){
            ArrayList given_sector = gridArr[sector];
            foreach (int index in given_sector){
                if (index!= particle_index){
                    float distance = magnitude(positions[index] - positions[particle_index]);

                    if (distance > 0.01f && distance <= smoothingRadius){
                        float2 direction = (positions[index] - positions[particle_index])/distance;
                        float deriv = smoothingKernelDerivativeNew(smoothingRadius, distance);
                        float sharedPressure = getSharedPressure(densities[index], densities[particle_index]);
                        pressureForce+= sharedPressure * direction * deriv * particle_mass / densities[index];
                    }
                    else{
                        float2 random_dir = new float2(Range(-1f, 1f), Range(-1f, 1f));
                        random_dir /= magnitude(random_dir);
                        float deriv = smoothingKernelDerivativeNew(smoothingRadius, distance);
                        float sharedPressure = getSharedPressure(densities[index], densities[particle_index]);
                        pressureForce += sharedPressure * deriv* random_dir * particle_mass / densities[index];

                    }
                }
            }
        }

        return pressureForce;
    }

    float getSharedPressure(float density1, float density2){
        float pressure_1 = densityToPressure(density1);
        float pressure_2 = densityToPressure(density2);
        return (pressure_1 + pressure_2)/2;
    }

    float densityToPressure(float density){
        float density_error = density - ideal_density;
        return density_error * pressure_factor;
    }

}
