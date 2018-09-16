# HybridRayTracing-Unity
Hybrid Ray Tracer made in Unity. 

![alt text](https://github.com/OpenRayTracing/HybridRayTracing-Unity/blob/master/Images/result16.png)

## Getting Started
Add HybridRayTracer.cs script to the main camera in the scene. The script has properties that can be changed from the Inspector. Use the example material to create metallic or glass materials. Unitys default material will be seen as diffuse material by the Hybrid Ray Tracer.

## Introduction 

### Per-Pixel linked list
A Per-Pixel linked list is used to store fragment information that are calculated from Unitys rendering pipeline. The infromation stored in each element is:
* color
* normal
* depth
* next (pointer to next element)
* v0 (Vertices of a triangle)
* v1
* v2
* mat (Material property)
* var (Variable to be used with metal or glass material)

Where v0/1/2 are obtained using a Geometry shader.

### Ray Tracing
The ray tracing step uses screen space ray tracing (from [kode80SSR](https://github.com/kode80/kode80SSR)) which visits pixels that are occupied by a given ray. Pixel position given by the screen space ray tracing is used to access the linked list for that pixel. Ray-Triangle intersection test is performed by iterating through the linked list using the vertices position stored in each element.
