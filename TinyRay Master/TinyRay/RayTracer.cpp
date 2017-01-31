/*---------------------------------------------------------------------
*
* Copyright © 2015  Minsi Chen
* E-mail: m.chen@derby.ac.uk
*
* The source is written for the Graphics I and II modules. You are free
* to use and extend the functionality. The code provided here is functional
* however the author does not guarantee its performance.
---------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <limits.h>


#if defined(WIN32) || defined(_WINDOWS)
#include <Windows.h>
#include <gl/GL.h>
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#endif

#include "RayTracer.h"
#include "Ray.h"
#include "Scene.h"
#include "Camera.h"
#include "perlin.h"


RayTracer::RayTracer()
{
	m_buffHeight = m_buffWidth = 0.0;
	m_renderCount = 0;
	SetTraceLevel(5);
	m_traceflag = (TraceFlags)(TRACE_AMBIENT | TRACE_DIFFUSE_AND_SPEC |
		TRACE_SHADOW | TRACE_REFLECTION | TRACE_REFRACTION);

}

RayTracer::RayTracer(int Width, int Height)
{
	m_buffWidth = Width;
	m_buffHeight = Height;
	m_renderCount = 0;
	SetTraceLevel(5);

	m_framebuffer = new Framebuffer(Width, Height);

	//default set default trace flag, i.e. no lighting, non-recursive
	m_traceflag = (TraceFlags)(TRACE_AMBIENT);
}

RayTracer::~RayTracer()
{
	delete m_framebuffer;
}

void RayTracer::DoRayTrace( Scene* pScene )
{
	Camera* cam = pScene->GetSceneCamera();
	
	Vector3 camRightVector = cam->GetRightVector();
	Vector3 camUpVector = cam->GetUpVector();
	Vector3 camViewVector = cam->GetViewVector();
	Vector3 centre = cam->GetViewCentre();
	Vector3 camPosition = cam->GetPosition();

	double sceneWidth = pScene->GetSceneWidth();
	double sceneHeight = pScene->GetSceneHeight();

	double pixelDX = sceneWidth / m_buffWidth;
	double pixelDY = sceneHeight / m_buffHeight;
	
	int total = m_buffHeight*m_buffWidth;
	int done_count = 0;
	
	Vector3 start;

	start[0] = centre[0] - ((sceneWidth * camRightVector[0])
		+ (sceneHeight * camUpVector[0])) / 2.0;
	start[1] = centre[1] - ((sceneWidth * camRightVector[1])
		+ (sceneHeight * camUpVector[1])) / 2.0;
	start[2] = centre[2] - ((sceneWidth * camRightVector[2])
		+ (sceneHeight * camUpVector[2])) / 2.0;
	
	Colour scenebg = pScene->GetBackgroundColour();

	if (m_renderCount == 0)
	{
		fprintf(stdout, "Trace start.\n");

		Colour colour;
//TinyRay on multiprocessors using OpenMP!!!
#pragma omp parallel for schedule (dynamic, 1) private(colour)
		for (int i = 0; i < m_buffHeight; i+=1) {
			for (int j = 0; j < m_buffWidth; j+=1) {

				//calculate the metric size of a pixel in the view plane (e.g. framebuffer)
				Vector3 pixel;

				pixel[0] = start[0] + (i + 0.5) * camUpVector[0] * pixelDY
					+ (j + 0.5) * camRightVector[0] * pixelDX;
				pixel[1] = start[1] + (i + 0.5) * camUpVector[1] * pixelDY
					+ (j + 0.5) * camRightVector[1] * pixelDX;
				pixel[2] = start[2] + (i + 0.5) * camUpVector[2] * pixelDY
					+ (j + 0.5) * camRightVector[2] * pixelDX;

				/*
				* setup first generation view ray
				* In perspective projection, each view ray originates from the eye (camera) position 
				* and pierces through a pixel in the view plane
				*/
				Ray viewray;
				viewray.SetRay(camPosition,	(pixel - camPosition).Normalise());
				
				double u = (double)j / (double)m_buffWidth;
				double v = (double)i / (double)m_buffHeight;

scenebg = pScene->GetBackgroundColour();

//trace the scene using the view ray
//default colour is the background colour, unless something is hit along the way
colour = this->TraceScene(pScene, viewray, scenebg, m_traceLevel);

/*
* Draw the pixel as a coloured rectangle
*/
m_framebuffer->WriteRGBToFramebuffer(colour, j, i);
			}
		}

		fprintf(stdout, "Done!!!\n");
		m_renderCount++;
	}
}

Colour RayTracer::TraceScene(Scene* pScene, Ray& ray, Colour incolour, int tracelevel, bool shadowray)
{
	RayHitResult result;

	Colour outcolour = incolour; //the output colour based on the ray-primitive intersection

	if (tracelevel <= 0) // Check for maximum depth of recursion
	{
		return outcolour;
	}

	std::vector<Light*> *light_list = pScene->GetLightList();
	Vector3 cameraPosition = pScene->GetSceneCamera()->GetPosition();

	//Intersect the ray with the scene
	//TODO: Scene::IntersectByRay needs to be implemented first
	result = pScene->IntersectByRay(ray);

	if (result.data) //the ray has hit something
	{
		//TODO:
		//1. Non-recursive ray tracing:
		//	 When a ray intersect with an objects, determine the colour at the intersection point
		//	 using CalculateLighting
		if (!shadowray) // Check if shadow ray and adjust surface colour accordingly
		{
			outcolour = CalculateLighting(light_list, &cameraPosition, &result);
		}
		else
		{
			if (!((Primitive*)result.data)->m_primtype == Primitive::PRIMTYPE_Plane)
				return (outcolour * 0.2);
		}

		//TODO: 2. The following conditions are for implementing recursive ray tracing
		if (m_traceflag & TRACE_REFLECTION)
		{
			//TODO: trace the reflection ray from the intersection point
			// Only reflect off primitives other than the planes
			if (!((Primitive*)result.data)->m_primtype == Primitive::PRIMTYPE_Plane)
			{
				Vector3 reflection_vec = ray.GetRay().Reflect(result.normal); // reflect about the normal to get the reflection vector
				Vector3 offset = reflection_vec * FLT_EPSILON * 50; // offset to reduce noise

				Ray reflection_ray;	// Create ray from the surface point for the reflected ray
				reflection_ray.SetRay(result.point + offset, reflection_vec);

				Colour reflectResult = TraceScene(pScene, reflection_ray, incolour, tracelevel-1, shadowray);
				outcolour = reflectResult * outcolour;
			}
		}

		if (m_traceflag & TRACE_REFRACTION)
		{
			//TODO: trace the refraction ray from the intersection point
			if (!((Primitive*)result.data)->m_primtype == Primitive::PRIMTYPE_Plane)
			{
				Vector3 refraction_vec = ray.GetRay().Refract(result.normal, 0.94);// 1.05);
				Vector3 offset = refraction_vec * FLT_EPSILON * 600; // offset to reduce noise

				Ray refraction_ray;
				refraction_ray.SetRay(result.point + offset, refraction_vec);

				Colour refractResult = TraceScene(pScene, refraction_ray, incolour, tracelevel-1, shadowray);
				outcolour = refractResult + outcolour;
			}
		}

		if (m_traceflag & TRACE_SHADOW)
		{
			//TODO: trace the shadow ray from the intersection point
			// Create new ray from intersection point to the light and check if it hits anything then set shadowray to true
			std::vector<Light*>::iterator lit_iter = light_list->begin();
			while (lit_iter != light_list->end()) // Iterate through all the lights in the scene
			{
				Vector3 ray = ((*lit_iter)->GetLightPosition() - result.point).Normalise(); // Vector from surface point towards light source
				Vector3 offset = ray * FLT_EPSILON * 10;

				Ray shadow_ray;
				shadow_ray.SetRay(result.point + offset, ray); //  Create ray from the surface point for the shadow ray

				outcolour = TraceScene(pScene, shadow_ray, outcolour, tracelevel-1, true);
				
				lit_iter++;
			}
		}
	}
		
	return outcolour;
}

Colour RayTracer::CalculateLighting(std::vector<Light*>* lights, Vector3* campos, RayHitResult* hitresult)
{
	Colour outcolour;
	std::vector<Light*>::iterator lit_iter = lights->begin();

	Primitive* prim = (Primitive*)hitresult->data;
	Material* mat = prim->GetMaterial();

	outcolour = mat->GetAmbientColour();
	
	//Generate the grid pattern on the plane
	if (((Primitive*)hitresult->data)->m_primtype == Primitive::PRIMTYPE_Plane)
	{
		int dx = hitresult->point[0]/2.0;
		int dy = hitresult->point[1]/2.0;
		int dz = hitresult->point[2]/2.0;

		if (dx % 2 || dy % 2 || dz % 2 )
		{
			outcolour = Vector3(0.1, 0.1, 0.1);
		}
		else
		{
			outcolour = mat->GetDiffuseColour();
		}
		return outcolour;
	}
	
	////Go through all lights in the scene
	////Note the default scene only has one light source
	if (m_traceflag & TRACE_DIFFUSE_AND_SPEC)
	{
		//TODO: Calculate and apply the lighting of the intersected object using the illumination model covered in the lecture
		//i.e. diffuse using Lambertian model, for specular, you can use either Phong or Blinn-Phong model
		while (lit_iter != lights->end()) // while iterator is not at its end
		{
			Colour surface_colour = mat->GetDiffuseColour();			// surface diffuse colour
			Colour light_colour = (*lit_iter)->GetLightColour();		// Intensity of light colour
			Vector3 normal = hitresult->normal;							// Get surface normal
			Vector3 light_position = (*lit_iter)->GetLightPosition();	// Position of the light source
			Vector3 surface_pt = hitresult->point;						// surface intersection point
			
			Vector3 light_vec = (light_position - surface_pt).Normalise(); // Calculate Light vector
			double theta = light_vec.DotProduct(normal);				// Calculate angle between light vector and the normal 

			// Diffuse Reflection - Lambertian Shading model (k*l)*cosθ 
			Colour resultDiffuse = (surface_colour * light_colour) * max(0, theta);

			Vector3 camera_vec = (*campos - surface_pt).Normalise();	// Camera direction vector
			
			// Specular Reflectance
			Vector3 half = (camera_vec + light_vec).Normalise();		// vector half way between light vector and the normal
			double phi1 = normal.DotProduct(half);						// angle between half vector and the normal
			Colour specular_colour = mat->GetSpecularColour();
			double specular_power = mat->GetSpecPower();

			// Blinn Phong model - (k*l)*(cos*φ')^n
			Colour resultSpecular = (specular_colour * light_colour) * pow(max(0, phi1), specular_power*2);

			Vector3 reflection_vec = light_vec.Reflect(normal);			// reflection of light vector in the normal
			double phi = -camera_vec.DotProduct(reflection_vec);		// angle between camera and reflection vector

			// Phong model - (k*l)*(cos*φ)^n
			resultSpecular = (specular_colour * light_colour) * pow(max(0, phi), specular_power); 

			// add them all together to get the colour of the surface point
			outcolour = resultDiffuse + resultSpecular + outcolour;

			lit_iter++;
		}
	}

	return outcolour;
}

