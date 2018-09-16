
Shader "Hidden/PerPixelLinkedList"
{
	Properties
	{
		_ID("ID", Int) = 0
		_Color("Color Tint", Color) = (1, 1, 1, 1)
		_MainTex("Main Tex", 2D) = "white" {}
		_NormalMap("Normal Map", 2D) = "white" {}
		//_BumpMap("Normal Map", 2D) = "bump" {}
		_Diffuse("Diffuse", Range(0,1)) = 0
		_Metal("Metal", Range(0,1)) = 0
		_Glass("Glass", Range(0,1)) = 0
		_Var("Var", Range(0,1)) = 0
		//_Layers("Layers", Int) = 2
	}
	SubShader
	{
		Tags{ "Queue" = "Transparent" }//"IgnoreProjector" = "True" "RenderType" = "Opaque" }
										//LOD 100

		Cull Off
		ZTest Always
		ZWrite Off




		Pass
		{
			//Tags{ "LightMode" = "ForwardBase" }

			/*
			Cull Off
			ZWrite Off
			ZTest Always

			*/



			CGPROGRAM
			#pragma target 5.0
			#pragma vertex vert
			#pragma geometry geo
			#pragma fragment frag

			#define MAX_FRAGMENTS 20

			#include "UnityCG.cginc"
			#include "Lighting.cginc"

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
				float3 normal : NORMAL;
				float4 tangent: TANGENT;
			};

			struct v2f
			{
				float4 vertex : SV_POSITION;
				float depth : DEPTH;
				float2 uv : TEXCOORD0;
				float4 screenPos: TEXCOORD1;

				nointerpolation float3 cameraPosition : TEXCOORD2;
				nointerpolation float3 v1 : TEXCOORD3;
				nointerpolation float3 v2 : TEXCOORD4;

				float3 normal: TEXCOORD5;
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


			float4x4 _V;
			float4x4 _P;
			float4x4 _VP;
			float4x4 _InvVP;

			float _nearClip;

			int _ID;
			float4 _Color;
			sampler2D _MainTex;
			float4 _MainTex_ST;
			sampler2D _BumpMap;
			int width;
			float4 resolution;
			float4 camPos;
			float4 camDir;
			float4 camUp;
			float4 camRight;

			int _Metal;
			int _Diffuse;
			int _Glass;
			float _Var;

			sampler2D _CameraDepthTexture;
			sampler2D _LastCameraDepthTexture;

			float4x4 _CameraInverseProjectionMatrix;

			RWStructuredBuffer<Node> list : register(u1);
			RWByteAddressBuffer listHead : register(u2);
			RWStructuredBuffer<int> counter : register(u3);


			v2f vert(appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.uv = TRANSFORM_TEX(v.uv, _MainTex);

				o.screenPos = ComputeScreenPos(o.vertex);
				o.depth = COMPUTE_DEPTH_01;
				o.normal = normalize(mul(UNITY_MATRIX_MV, v.normal));

				o.cameraPosition = mul(UNITY_MATRIX_MV, v.vertex).xyz;

				return o;
			}

			[maxvertexcount(3)]
			void geo(triangle v2f input[3], inout TriangleStream<v2f> OutputStream)
			{
				v2f o = (v2f)0;

				o.normal = input[0].normal;
				o.vertex = input[0].vertex;
				o.uv = input[0].uv;
				o.depth = input[0].depth;
				o.screenPos = input[0].screenPos;
				o.cameraPosition = input[0].cameraPosition;
				o.v1 = input[1].cameraPosition;
				o.v2 = input[2].cameraPosition;
				OutputStream.Append(o);

				o.normal = input[1].normal;
				o.vertex = input[1].vertex;
				o.uv = input[1].uv;
				o.depth = input[1].depth;
				o.screenPos = input[1].screenPos;
				o.cameraPosition = input[1].cameraPosition;
				o.v1 = input[2].cameraPosition;
				o.v2 = input[0].cameraPosition;
				OutputStream.Append(o);

				o.normal = input[2].normal;
				o.vertex = input[2].vertex;
				o.uv = input[2].uv;
				o.depth = input[2].depth;
				o.screenPos = input[2].screenPos;
				o.cameraPosition = input[2].cameraPosition;
				o.v1 = input[0].cameraPosition;
				o.v2 = input[1].cameraPosition;
				OutputStream.Append(o);


				OutputStream.RestartStrip();
			}

			void frag(v2f i)
			{
				float4 texColor = tex2D(_MainTex, i.uv);

				float4 col = texColor * _Color;

				uint newFragmentAddress = counter.IncrementCounter();

				float2 uv = i.screenPos.xy / i.screenPos.w;
				uint2 screenpos = 0.5 * (uv + 1.0) *_ScreenParams.xy;

				uint bufferPos = 4 * ((width*screenpos.y) + screenpos.x);
				uint oldFragmentAddress;

				listHead.InterlockedExchange(bufferPos, newFragmentAddress, oldFragmentAddress);

				Node n;
				n.color = col.rgb;
				n.normal = normalize(i.normal);
				n.depth = i.depth;
				n.next = oldFragmentAddress;
				n.v0 = i.cameraPosition;
				n.v1 = i.v1;
				n.v2 = i.v2;
				n.mat = float3(_Diffuse, _Metal, _Glass);
				n.var = _Var;

				list[newFragmentAddress] = n;

				return;


				/* TEST SORTING
				if (oldFragmentAddress != 0xffffffff) {
				uint currentNode = oldFragmentAddress;
				Node tempList[MAX_FRAGMENTS];
				uint num = 0;
				while (currentNode != 0xffffffff && num < MAX_FRAGMENTS) {
				tempList[num] = list[currentNode];
				currentNode = tempList[num].next;
				num++;
				}

				for (uint j = 0; j < num; j++) {
				if (n.depth < tempList[j].depth) {
				//update headlist, med load gjør vi det bare hvis new e nærmest
				if (j != 0) {
				//den før og etter hvis etter finnes
				uint after = tempList[j].mem;//list[tempList[j - 1].mem].next;
				list[tempList[j - 1].mem].next = n.mem;
				if (j < MAX_FRAGMENTS) {//num-1){
				n.next = after;
				}
				else {
				n.next = 0xffffffff;
				}
				}
				else {
				listHead.InterlockedExchange(bufferPos, newFragmentAddress, old);
				n.next = old;
				}

				list[newFragmentAddress] = n;
				break;
				}
				else {
				if (j == num - 1 && j < MAX_FRAGMENTS) {
				n.next = 0xffffffff;
				linkedList[tempList[j].mem].next = n.mem;
				linkedList[newFragmentAddress] = n;
				break;
				}
				}
				}
				}
				else {
				listHead.InterlockedExchange(bufferPos, newFragmentAddress, old);
				n.next = old;
				linkedList[newFragmentAddress] = n;
				}
				return
				*/
			}
			ENDCG
		}
	}
}
