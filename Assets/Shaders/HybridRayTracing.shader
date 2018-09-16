// Upgrade NOTE: replaced '_Object2World' with 'unity_ObjectToWorld'
// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

Shader "Hidden/HybridRayTracing"
{
	Properties
	{
		_Lambertian("Lambertian", Range(0,1)) = 1
		_Dielectric("Dielectric", Range(0,1)) = 0
		_Metal("Metal", Range(0,1)) = 0
		_MainTex ("Texture", 2D) = "white" {}
		_PrevFrame("Texture", 2D) = "white" {}
		_Color("Color", Color) = (1,1,1,1)
		_Cube("Reflection Map", Cube) = "" {}
	}
	SubShader
	{
		Tags { "RenderType" = "Opaque" }//"RenderType" = "Opaque" "Queue" = "Transparent"}
		//ZTest off Cull Off ZWrite Off
		//Blend Off
		//Cull Off
		//ZTest Always
		//ZWrite Off

		Pass
		{
			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma fragment frag
			
			#include "UnityCG.cginc"

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
			};

			struct v2f
			{
				float2 uv : TEXCOORD0;
				float4 vertex : SV_POSITION;
				float4 screenPos: TEXCOORD1;
				float4 cameraRay : TEXCOORD2;
			};

			struct Node {
				float3 color;
				float3 normal;
				float depth;
				uint next;
				float3 v0;
				float3 v1;
				float3 v2;
				float3 mat;
				float var;
			};

			struct Ray {
				float3 origin;
				float3 direction;
				void make(float3 orig, float3 dir) {
					origin = orig;
					direction = dir;
				}
			};

			sampler2D _MainTex;
			float4 _MainTex_ST;

			sampler2D _CameraGBufferTexture0;	// Diffuse color (RGB), unused (A)
			sampler2D _CameraGBufferTexture1;	// Specular color (RGB), roughness (A)
			sampler2D _CameraGBufferTexture2;	// World space normal (RGB), unused (A)
			sampler2D _CameraGBufferTexture3;	// ARGBHalf (HDR) format: Emission + lighting + lightmaps + reflection probes buffer

			sampler2D _PrevFrame;

			int _Refinement;

			samplerCUBE _Cube;

			int _Lambertian;
			int _Dielectric;
			int _Metal;

			float4 _Color;

			StructuredBuffer<Node> list : register(t1);
			Buffer<uint> head : register(t2);

			float4x4 _V;
			float4x4 _P;
			float4x4 _VP;
			float4x4 _InvVP;
			float4x4 _InvP;
			float4x4 _InvV;

			float _nearClip;

			float _RayLength;
			float _RayThickness;
			int _RayNum;
			int _RayBounce;

			int width;
			float4 resolution;
			float4 camPos;
			float4 camDir;
			float4 camUp;
			float4 camRight;


			static float epsilon = 1e-8;

			sampler2D _CameraDepthTexture;
			sampler2D _BackFaceDepthTex;

			int _FrameCount;
			int _RealFrameCount;
			int _Raytrace;
			int _SSR;
			int _RPP;

			float4x4 _CameraProjectionMatrix;			// projection matrix that maps to screen pixels (not NDC)
			float4x4 _CameraInverseProjectionMatrix;	// inverse projection matrix (NDC to camera space)
			float _Iterations;							// maximum ray iterations
			float _BinarySearchIterations;				// maximum binary search refinement iterations
			float _PixelZSize;							// Z size in camera space of a pixel in the depth buffer
			float _PixelStride;							// number of pixels per ray step close to camera
			float _PixelStrideZCuttoff;					// ray origin Z at this distance will have a pixel stride of 1.0
			float _MaxRayDistance;						// maximum distance of a ray
			float _ScreenEdgeFadeStart;					// distance to screen edge that ray hits will start to fade (0.0 -> 1.0)
			float _EyeFadeStart;						// ray direction's Z that ray hits will start to fade (0.0 -> 1.0)
			float _EyeFadeEnd;							// ray direction's Z that ray hits will be cut (0.0 -> 1.0)

			float4x4 _NormalMatrix;
			float2 _RenderBufferSize;
			float2 _OneDividedByRenderBufferSize;		// Optimization: removes 2 divisions every itteration


			v2f vert(appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.uv = TRANSFORM_TEX(v.uv, _MainTex);
				o.screenPos = ComputeScreenPos(o.vertex);
				//o.rayDir = -WorldSpaceViewDir(v.vertex);
				float4 cameraRay = float4(o.uv * 2.0 - 1.0, 1.0, 1.0);
				cameraRay = mul(_CameraInverseProjectionMatrix, cameraRay);
				o.cameraRay = cameraRay / cameraRay.w;
				return o;
			}

			float rand_1_05(in float2 uv)
			{
				float2 noise = (frac(sin(dot(uv, float2(12.9898, 78.233)*2.0)) * 43758.5453));
				return abs(noise.x + noise.y) * 0.5;
			}

			float noise(float2 seed)
			{
				return frac(sin(dot(seed.xy, float2(12.9898, 78.233))) * 43758.5453);
			}

			float3 random_in_unit_sphere(float2 uv) {
				float3 p;
				uint i = 1;
				do {
					uint x = i + 1;
					uint y = i * (i + 2);
					uint z = (i*(i + 2)) / 2;
					i++;
					p = 2.0*float3(rand_1_05(uv + x), rand_1_05(uv + y), rand_1_05(uv + z)) - float3(1, 1, 1);
				} while (sqrt(length(p)) >= 1.0);
				return p;
			}

			float3 random_in_unit_disk(float2 uv) {
				float3 p;
				int i = 1;
				do {
					int x = i + 1;
					int y = i * (i + 2);
					i++;
					p = 2.0*float3(rand_1_05(x + uv), rand_1_05(y + uv), 0) - float3(1, 1, 0);
				} while (dot(p, p) >= 1.0);
				return p;
			}

			class camera {
				float3 lower_left_corner;
				float3 horizontal;
				float3 vertical;
				float3 origin;
				float3 u, v, w;
				float lens_radius;


				void make(float3 lookfrom, float3 lookat, float3 vup, float vfov, float aspect, float aperture, float focus_dist) {
					lens_radius = aperture / 2;
					float theta = vfov * 3.14159265 / 180;
					float half_height = tan(theta / 2);
					float half_width = aspect * half_height;
					origin = lookfrom;
					w = normalize(lookfrom - lookat);
					u = normalize(cross(vup, w));
					v = cross(w, u);
					lower_left_corner = origin - half_width * focus_dist*u - half_height * focus_dist*v - focus_dist * w;
					horizontal = 2 * half_width*focus_dist*u;
					vertical = 2 * half_height*focus_dist*v;
				}

				Ray getRay(float s, float t) {
					float3 rd = lens_radius * random_in_unit_disk(float2(s, t));
					float3 offset = u * rd.x + v * rd.y;
					Ray r;
					r.make(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset);
					return r;
				}
			};
			
			bool rayTriangleIntersection(float3 origin, float3 direction, float3 v0, float3 v1, float3 v2, out float t) {

				float3 e1 = v1 - v0;
				float3 e2 = v2 - v0;

				float3 N = cross(e1, e2);
				float area2 = length(N);

				float NdotRayDirection = dot(N, direction);
				if (abs(NdotRayDirection) < epsilon) return false;


				float d = dot(N, v0);
				t = (d - dot(N, origin)) / NdotRayDirection;
				if (t < 0) return false; 

				float3 P = origin + t * direction;


				float3 C;
				float3 edge0 = v1 - v0;
				float3 vp0 = P - v0;
				C = cross(edge0, vp0);
				if (dot(N, C) < 0) return false;

				float3 edge1 = v2 - v1;
				float3 vp1 = P - v1;
				C = cross(edge1, vp1);
				if (dot(N, C) < 0) return false; 

				float3 edge2 = v0 - v2;
				float3 vp2 = P - v2;
				C = cross(edge2, vp2);
				if (dot(N, C) < 0) return false;

				return true;
			}

			float ComputeDepth(float4 clippos)
			{
				/*
				#if defined(SHADER_TARGET_GLSL) || defined(SHADER_API_GLES) || defined(SHADER_API_GLES3)
				return (clippos.z / clippos.w) * 0.5 + 0.5;
				#else
				return clippos.z / clippos.w;
				#endif
				*/
				return clippos.z / clippos.w;
			}

			bool raytraceTriangles(in uint index, in Ray r, inout float mint, inout float maxt, inout Node currentNode) {
				
				float t;
				bool found = false;

				while (index != 0xffffffff) {
					Node node = list[index];
					float3 v0, v1, v2;
					v0 = node.v0;
					v1 = node.v1;
					v2 = node.v2;

					bool intersec = rayTriangleIntersection(r.origin, r.direction, v0, v1, v2, t);

					if (intersec) {
						if (t < maxt && t > mint) {
							maxt = t;
							found = true;
							currentNode = node;
						}
					}

					index = node.next;
				}

				return found;
			}
			
			float4 doTracing() {

			}

			inline bool rayIntersectsDepthBF2(float zA, float zB, float2 uv)
			{
				uint2 screenpos = (0.5 * (uv + 1.0)) * _ScreenParams.xy;
				uint bufferPos = width * screenpos.y + screenpos.x;
				uint index = head.Load(bufferPos);
				return index != -1;
				//Node n = list[index];
				//return  n.depth < 1;
			}
			
			
			//-------------------------------------------------------------------------------------------------------------------------------------
			//  Copyright (c) 2015, Ben Hopkins (kode80)
			//  All rights reserved.
			//  
			//  Redistribution and use in source and binary forms, with or without modification, 
			//  are permitted provided that the following conditions are met:
			//  
			//  1. Redistributions of source code must retain the above copyright notice, 
			//     this list of conditions and the following disclaimer.
			//  
			//  2. Redistributions in binary form must reproduce the above copyright notice, 
			//     this list of conditions and the following disclaimer in the documentation 
			//     and/or other materials provided with the distribution.
			//  
			//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
			//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
			//  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
			//  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
			//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
			//  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
			//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
			//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
			//  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
			inline float3 ScreenSpaceToViewSpace(float3 cameraRay, float depth)
			{
				return (cameraRay * depth);
			}

			inline void swapIfBigger(inout float aa, inout float bb)
			{
				if (aa > bb)
				{
					float tmp = aa;
					aa = bb;
					bb = tmp;
				}
			}

			inline bool rayIntersectsDepthBF(float zA, float zB, float2 uv)
			{
				float4 uv4 = float4(uv, 0.0, 0.0);
				float cameraZ = Linear01Depth(tex2Dlod(_CameraDepthTexture, uv4).r) * -_ProjectionParams.z;
				float backZ = tex2Dlod(_BackFaceDepthTex, uv4).r * -_ProjectionParams.z;

				return zB <= cameraZ && zA >= backZ - _PixelZSize;
			}

			inline float distanceSquared(float2 a, float2 b) { a -= b; return dot(a, a); }

			// Trace a ray in screenspace from rayOrigin (in camera space) pointing in rayDirection (in camera space)
			// using jitter to offset the ray based on (jitter * _PixelStride).
			//
			// Returns true if the ray hits a pixel in the depth buffer
			// and outputs the hitPixel (in UV space), the hitPoint (in camera space) and the number
			// of iterations it took to get there.
			//
			// Based on Morgan McGuire & Mike Mara's GLSL implementation:
			// http://casual-effects.blogspot.com/2014/08/screen-space-ray-tracing.html
			inline bool traceScreenSpaceRay(
				in float3 rayOrigin,
				in float3 rayDirection,
				float jitter,
				out float2 hitPixel,
				out float3 hitPoint,
				out float iterationCount, 
				bool debugHalf, 
				inout Ray ray,
				in float mint,
				inout float maxt,
				inout Node n,
				bool checkAll
				)
			{
				// Clip to the near plane    
				float rayLength = ((rayOrigin.z + rayDirection.z * _MaxRayDistance) > -_ProjectionParams.y) ?
					(-_ProjectionParams.y - rayOrigin.z) / rayDirection.z : _MaxRayDistance;
				float3 rayEnd = rayOrigin + rayDirection * rayLength;

				// Project into homogeneous clip space
				float4 H0 = mul(_CameraProjectionMatrix, float4(rayOrigin, 1.0));
				float4 H1 = mul(_CameraProjectionMatrix, float4(rayEnd, 1.0));

				float k0 = 1.0 / H0.w, k1 = 1.0 / H1.w;

				// The interpolated homogeneous version of the camera-space points  
				float3 Q0 = rayOrigin * k0, Q1 = rayEnd * k1;

				// Screen-space endpoints
				float2 P0 = H0.xy * k0, P1 = H1.xy * k1;

				// If the line is degenerate, make it cover at least one pixel
				// to avoid handling zero-pixel extent as a special case later
				P1 += (distanceSquared(P0, P1) < 0.0001) ? 0.01 : 0.0;

				float2 delta = P1 - P0;

				// Permute so that the primary iteration is in x to collapse
				// all quadrant-specific DDA cases later
				bool permute = false;
				if (abs(delta.x) < abs(delta.y)) {
					// This is a more-vertical line
					permute = true; delta = delta.yx; P0 = P0.yx; P1 = P1.yx;
				}

				float stepDir = sign(delta.x);
				float invdx = stepDir / delta.x;

				// Track the derivatives of Q and k
				float3  dQ = (Q1 - Q0) * invdx;
				float dk = (k1 - k0) * invdx;
				float2  dP = float2(stepDir, delta.y * invdx);

				// Calculate pixel stride based on distance of ray origin from camera.
				// Since perspective means distant objects will be smaller in screen space
				// we can use this to have higher quality reflections for far away objects
				// while still using a large pixel stride for near objects (and increase performance)
				// this also helps mitigate artifacts on distant reflections when we use a large
				// pixel stride.
				float strideScaler = 1.0 - min(1.0, -rayOrigin.z / _PixelStrideZCuttoff);
				float pixelStride = _PixelStride;//1.0 + strideScaler * _PixelStride;

				// Scale derivatives by the desired pixel stride and then
				// offset the starting values by the jitter fraction
				dP *= pixelStride; dQ *= pixelStride; dk *= pixelStride;
				P0 += dP * jitter; Q0 += dQ * jitter; k0 += dk * jitter;

				float i, zA = 0.0, zB = 0.0;

				// Track ray step and derivatives in a float4 to parallelize
				float4 pqk = float4(P0, Q0.z, k0);
				float4 dPQK = float4(dP, dQ.z, dk);
				bool intersect = false;

				for (i = 0; i<_Iterations && intersect == false; i++)
				{
					pqk += dPQK;

					zA = zB;
					zB = (dPQK.z * 0.5 + pqk.z) / (dPQK.w * 0.5 + pqk.w);
					swapIfBigger(zB, zA);

					hitPixel = permute ? pqk.yx : pqk.xy;
					hitPixel *= _OneDividedByRenderBufferSize;
					
					if (checkAll){//_Raytrace) {
						bool inter = rayIntersectsDepthBF2(zA, zB, hitPixel);

						if (inter) {
							uint2 screenpos = (0.5 * (hitPixel + 1.0)) * _ScreenParams.xy;
							uint bufferPos = width * screenpos.y + screenpos.x;
							uint index = head.Load(bufferPos);
							//n = list[index];
							//maxt = 1000;
							//ray = { rayOrigin, rayDirection };
							intersect = raytraceTriangles(index, ray, mint, maxt, n);
						}
					}
					else {
						intersect = rayIntersectsDepthBF(zA, zB, hitPixel);
					}
					
					
					//intersect = rayIntersectsDepthBF(zA, zB, hitPixel);
				}
				
				// Binary search refinement
				/*
				if (pixelStride > 1.0 && intersect)
				{
					pqk -= dPQK;
					dPQK /= pixelStride;

					float originalStride = pixelStride * 0.5;
					float stride = originalStride;

					zA = pqk.z / pqk.w;
					zB = zA;

					for (float j = 0; j<_BinarySearchIterations; j++)
					{
						pqk += dPQK * stride;

						zA = zB;
						zB = (dPQK.z * -0.5 + pqk.z) / (dPQK.w * -0.5 + pqk.w);
						swapIfBigger(zB, zA);

						hitPixel = permute ? pqk.yx : pqk.xy;
						hitPixel *= _OneDividedByRenderBufferSize;

						originalStride *= 0.5;
						stride = rayIntersectsDepthBF(zA, zB, hitPixel) ? -originalStride : originalStride;
					}
				}
				*/

				Q0.xy += dQ.xy * i;
				Q0.z = pqk.z;
				hitPoint = Q0 / pqk.w;
				iterationCount = i;

				return intersect;
			}
			//-------------------------------------------------------------------------------------------------------------------------------------

			bool refra(float3 v, float3 n, float ni_over_nt, out float3 refracted) {
				float3 uv = normalize(v);
				float dt = dot(uv, n);
				float discriminant = 1.0 - ni_over_nt * ni_over_nt*(1 - dt * dt);
				if (discriminant > 0) {
					refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
					return true;
				}
				else {
					return false;
				}
			}

			float3 refractTest(in float3 I, in float3 N, in float ior) {
				float cosi = clamp(-1, 1, dot(I, N));
				float etai = 1, etat = ior;
				float3 n = N;
				if (cosi < 0) {
					cosi = -cosi; 
				}
				else { 
					float temp = etai;
					etai = etat;
					etat = temp;
					n = -N; 
				}
				float eta = etai / etat;
				float k = 1 - eta * eta * (1 - cosi * cosi);
				return k < 0 ? 0 : eta * I + (eta * cosi - sqrt(k)) * n;
			}

			float schlick(float cosine, float ref_idx) {
				float r0 = (1 - ref_idx) / (1 + ref_idx);
				r0 *= r0;
				return r0 + (1 - r0)*pow((1 - cosine), 5);
			}

			void fresnel(in float3 I, in float3 N, in float ior, inout float kr) {
				float cosi = clamp(-1, 1, dot(I, N));
				float etai = 1, etat = ior;
				if (cosi > 0) {
					float temp = etai;
					etai = etat;
					etat = temp;
				}

				float sint = etai / etat * sqrt(max(0.f, 1 - cosi * cosi));

				if (sint >= 1) {
					kr = 1;
				}
				else {
					float cost = sqrt(max(0.f, 1 - sint * sint));
					cosi = abs(cosi);
					float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
					float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
					kr = (Rs * Rs + Rp * Rp) / 2;
				}
			}

			float3 reflect(in float3 I, in float3 N)
			{
				return I - 2.0 * dot(I, N) * N;
			}

			float3 color(inout Ray r, float2 uv, float2 uuvv) {
				float3 col = float3(1, 1, 1);

				uint2 screenpos = (0.5 * (uv + 1.0)) * _ScreenParams.xy;
				uint bufferPos = width * screenpos.y + screenpos.x;
				uint index = head.Load(bufferPos);

				Node n = list[index];
				float maxt = 1000;
				float mint = 0.001;
				

				bool hit = raytraceTriangles(index, r, mint, maxt, n);
				
				if (!hit) {
					return float3(tex2D(_MainTex, uv).xyz);
				}
				
				float kr;

				float2 hitPixel;
				float3 hitPoint;
				int iterationCount;
				float2 uv2 = uv * _RenderBufferSize;
				float c = (uv2.x + uv2.y) * 0.25;
				float jitter = fmod(c, 1.0);
				int maxbounce = _RayBounce;
				while (hit && maxbounce > 0) {

					maxbounce--;
					hit = false;

					col *= n.color;
					

					float3 normal = n.normal;
					
					r.origin = r.origin + r.direction*maxt;
					if (n.mat.y == 1) {
						float3 reflected = reflect(normalize(r.direction), normal); 
						float fuzz = n.var;
						r.direction = reflected + fuzz * random_in_unit_sphere(uuvv);
						maxt = 1000;
						hit = traceScreenSpaceRay(r.origin, r.direction, jitter, hitPixel, hitPoint, iterationCount, uv.x > 0.5, r, mint, maxt, n, true);
					}
					else if (n.mat.z == 1) {
						float ref_idx = n.var;
						float3 outward_normal;
						float3 reflected = reflect(normalize(r.direction), normal);
						float ni_over_nt;
						float3 refracted;
						float reflect_prob;
						float cosine;
						if (dot(r.direction, normal) > 0) {
							outward_normal = -normal;
							ni_over_nt = ref_idx;
							cosine = ref_idx * dot(r.direction, normal) / length(r.direction);
						}
						else {
							outward_normal = normal;
							ni_over_nt = 1.0 / ref_idx;
							cosine = -dot(r.direction, normal) / length(r.direction);
						}

						if (refra(r.direction, outward_normal, ni_over_nt, refracted)) {
							reflect_prob = schlick(cosine, ref_idx);
						}
						else {
							r.direction = reflected;
							reflect_prob = 1.0;
						}
						if (rand_1_05(uuvv) < reflect_prob) {
							r.direction = reflected;
						}
						else {
							r.direction = refracted;
						}
						
						maxt = 1000;
						hit = traceScreenSpaceRay(r.origin, r.direction, jitter, hitPixel, hitPoint, iterationCount, uv.x > 0.5, r, mint, maxt, n, true);
						
					}
					else {
						float3 target = r.origin + normal + random_in_unit_sphere(uuvv);
						r.direction = target - r.origin;

						maxt = 1000;
						hit = traceScreenSpaceRay(r.origin, r.direction, jitter, hitPixel, hitPoint, iterationCount, uv.x > 0.5, r, mint, maxt, n, true);

						/*ONLY USE DEPTH TEST AND TEST ONCE
						maxt = 1000;
						hit = traceScreenSpaceRay(r.origin, r.direction, jitter, hitPixel, hitPoint, iterationCount, uv.x > 0.5, r, mint, maxt, n, true);
						if (hit) {
							screenpos = (0.5 * (hitPixel + 1.0)) * _ScreenParams.xy;
							bufferPos = width * screenpos.y + screenpos.x;
							index = head.Load(bufferPos);
							n = list[index];
							maxt = 1000;
							hit = raytraceTriangles(index, r, mint, maxt, n);
						}
						*/
					}
					
					
				}

				if (hit && maxbounce == 0) {
					col = float3(0, 0, 0);
				}
				float3 dd = mul(_InvV, float4(r.direction.xyz, 0)).xyz;
				float4 skyData = UNITY_SAMPLE_TEXCUBE(unity_SpecCube0, dd);
				float3 skyColor = DecodeHDR(skyData, unity_SpecCube0_HDR);


				return col * skyColor;
			}
			
			float4 frag (v2f i) : SV_Target
			{
				float2 uv = i.screenPos.xy / i.screenPos.w;

				
				float3 camRay = normalize(i.cameraRay.xyz);
				
				
				float3 lookAt = camRay * 10;
				float dist_to_focus = length(float3(0, 0, 0) - lookAt);
				camera cam;
				cam.make(float3(0, 0, 0), lookAt, float3(0, 1, 0), 0, 0, 0, dist_to_focus);

				float seed1 = uv.x + _Time.y;
				float seed2 = uv.y + _Time.y;
				float3 col = float3(0, 0, 0);
				for (int j = 0; j < _RPP; j++) {
					float uu = (uv.x + rand_1_05(seed1)); //(camRay.x + rand_1_05(seed1 + j) / 1000.0f);//(i.uv.x + rand_1_05(seed1 + j) / 400.0f);
					float vv = (uv.y + rand_1_05(seed2)); //(camRay.y + rand_1_05(seed2 + j * j) / 1000.0f);//(i.uv.y + rand_1_05(seed2 + j * j) / 400.0f);
																			//return float4(uu, vv, 0, 1);
					Ray r = cam.getRay(uu, vv); //{ float3(0,0,0), camRay.xyz };
					col += color(r, float2(uv), float2(uu, vv));
					
					
				}

				if (_FrameCount > _RPP && _Refinement == 1) {

					int time = _FrameCount;
					float4 prevCol = tex2D(_PrevFrame, i.uv);
					float3 avgColor = (col + prevCol.xyz*time) / (time + _RPP);
					col = avgColor;
				}
				else {
					col /= float(_RPP);
				}

				return float4(col, 1);
				
				
			}
			ENDCG
		}
	}
}
