// URP lit terrain shader that colors the surface by world-space elevation and slope.
// Five height bands (ocean → beach → grass → rock → snow), each tunable in the material
// Inspector, with smooth transitions controlled by a shared blend width. A separate slope
// override pushes steep faces toward the rock color regardless of altitude, so cliffs and
// mountainsides read correctly even at low elevations.
//
// Uses the world-space Y of each fragment, which means the same material works for every
// chunk regardless of where its GameObject sits (per-chunk meshes are in world coords; the
// shared flat-ocean tile is at WorldCenter — both end up with the same Y at sea level).

Shader "Custom/TerrainElevation"
{
    Properties
    {
        [Header(Height Bands)]
        _ColorOcean ("Ocean Color",      Color) = (0.08, 0.18, 0.35, 1)
        _ColorBeach ("Beach Color",      Color) = (0.85, 0.78, 0.55, 1)
        _ColorGrass ("Grass Color",      Color) = (0.30, 0.55, 0.25, 1)
        _ColorRock  ("Rock Color",       Color) = (0.45, 0.40, 0.38, 1)
        _ColorSnow  ("Snow Color",       Color) = (0.95, 0.95, 0.98, 1)

        [Header(Height Thresholds in world units)]
        _BeachHeight ("Beach Start",  Float) = 0.2
        _GrassHeight ("Grass Start",  Float) = 1.5
        _RockHeight  ("Rock Start",   Float) = 8.0
        _SnowHeight  ("Snow Start",  Float) = 15.0
        _BandBlend   ("Band Blend Width", Float) = 0.6

        [Header(Slope Override)]
        [Toggle] _SlopeEnable    ("Enable Slope Override", Float) = 1
        _SlopeThreshold ("Slope Threshold (0=flat, 1=vertical)", Range(0,1)) = 0.45
        _SlopeBlend     ("Slope Blend Width",                    Range(0,1)) = 0.15

        [Header(Lighting)]
        _AmbientBoost ("Ambient Boost", Range(0,1)) = 0.25
    }

    SubShader
    {
        Tags
        {
            "RenderType"        = "Opaque"
            "RenderPipeline"    = "UniversalPipeline"
            "Queue"             = "Geometry"
            "IgnoreProjector"   = "True"
        }
        LOD 100

        // ------------- Forward lit pass -------------
        Pass
        {
            Name "ForwardLit"
            Tags { "LightMode" = "UniversalForward" }

            HLSLPROGRAM
            #pragma vertex   Vert
            #pragma fragment Frag
            #pragma shader_feature_local _SLOPEENABLE_ON
            #pragma multi_compile _ _MAIN_LIGHT_SHADOWS _MAIN_LIGHT_SHADOWS_CASCADE _MAIN_LIGHT_SHADOWS_SCREEN
            #pragma multi_compile_fragment _ _SHADOWS_SOFT
            #pragma multi_compile_fog

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Shadows.hlsl"

            CBUFFER_START(UnityPerMaterial)
                float4 _ColorOcean;
                float4 _ColorBeach;
                float4 _ColorGrass;
                float4 _ColorRock;
                float4 _ColorSnow;
                float  _BeachHeight;
                float  _GrassHeight;
                float  _RockHeight;
                float  _SnowHeight;
                float  _BandBlend;
                float  _SlopeEnable;
                float  _SlopeThreshold;
                float  _SlopeBlend;
                float  _AmbientBoost;
            CBUFFER_END

            struct Attributes
            {
                float4 positionOS : POSITION;
                float3 normalOS   : NORMAL;
            };

            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 positionWS : TEXCOORD0;
                float3 normalWS   : TEXCOORD1;
                float  fogFactor  : TEXCOORD2;
            };

            Varyings Vert(Attributes IN)
            {
                Varyings OUT;
                VertexPositionInputs pos = GetVertexPositionInputs(IN.positionOS.xyz);
                VertexNormalInputs   nrm = GetVertexNormalInputs(IN.normalOS);
                OUT.positionCS = pos.positionCS;
                OUT.positionWS = pos.positionWS;
                OUT.normalWS   = nrm.normalWS;
                OUT.fogFactor  = ComputeFogFactor(pos.positionCS.z);
                return OUT;
            }

            // Smoothly blend `c` toward `next` as `h` crosses `start` (over a `blend` window).
            float3 BlendBand(float3 c, float3 next, float h, float start, float blend)
            {
                float t = smoothstep(start - blend, start + blend, h);
                return lerp(c, next, t);
            }

            half4 Frag(Varyings IN) : SV_Target
            {
                float3 N = normalize(IN.normalWS);
                float  h = IN.positionWS.y;

                // 1) Height-banded base color.
                float3 col = _ColorOcean.rgb;
                col = BlendBand(col, _ColorBeach.rgb, h, _BeachHeight, _BandBlend);
                col = BlendBand(col, _ColorGrass.rgb, h, _GrassHeight, _BandBlend);
                col = BlendBand(col, _ColorRock.rgb,  h, _RockHeight,  _BandBlend);
                col = BlendBand(col, _ColorSnow.rgb,  h, _SnowHeight,  _BandBlend);

                // 2) Slope override — steep faces become rock-colored. slope = 0 on horizontal
                //    ground, 1 on a vertical wall.
                #ifdef _SLOPEENABLE_ON
                {
                    float slope = 1.0 - saturate(N.y);
                    float t = smoothstep(_SlopeThreshold - _SlopeBlend,
                                         _SlopeThreshold + _SlopeBlend, slope);
                    col = lerp(col, _ColorRock.rgb, t);
                }
                #endif

                // 3) Standard URP lighting (main directional + soft ambient SH).
                Light mainLight  = GetMainLight(TransformWorldToShadowCoord(IN.positionWS));
                float NdotL      = saturate(dot(N, mainLight.direction));
                float3 direct    = col * mainLight.color * NdotL * mainLight.shadowAttenuation;
                float3 ambient   = col * (SampleSH(N) + _AmbientBoost);
                float3 lit       = direct + ambient;

                lit = MixFog(lit, IN.fogFactor);
                return half4(lit, 1.0);
            }
            ENDHLSL
        }

        // ------------- Shadow caster pass -------------
        // So terrain casts shadows onto itself and other receivers (boat, NPCs).
        Pass
        {
            Name "ShadowCaster"
            Tags { "LightMode" = "ShadowCaster" }

            ZWrite On
            ZTest LEqual
            ColorMask 0

            HLSLPROGRAM
            #pragma vertex   ShadowVert
            #pragma fragment ShadowFrag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Shadows.hlsl"

            float3 _LightDirection;
            float3 _LightPosition;

            struct AttributesS
            {
                float4 positionOS : POSITION;
                float3 normalOS   : NORMAL;
            };

            struct VaryingsS
            {
                float4 positionCS : SV_POSITION;
            };

            VaryingsS ShadowVert(AttributesS IN)
            {
                VaryingsS OUT;
                float3 positionWS = TransformObjectToWorld(IN.positionOS.xyz);
                float3 normalWS   = TransformObjectToWorldNormal(IN.normalOS);
                float4 positionCS = TransformWorldToHClip(ApplyShadowBias(positionWS, normalWS, _LightDirection));

                #if UNITY_REVERSED_Z
                    positionCS.z = min(positionCS.z, UNITY_NEAR_CLIP_VALUE);
                #else
                    positionCS.z = max(positionCS.z, UNITY_NEAR_CLIP_VALUE);
                #endif

                OUT.positionCS = positionCS;
                return OUT;
            }

            half4 ShadowFrag(VaryingsS IN) : SV_Target { return 0; }
            ENDHLSL
        }

        // ------------- DepthOnly pass -------------
        // Needed for SSAO, depth texture, and the project's DepthNormalsFeature.
        Pass
        {
            Name "DepthOnly"
            Tags { "LightMode" = "DepthOnly" }

            ZWrite On
            ColorMask 0

            HLSLPROGRAM
            #pragma vertex   DepthVert
            #pragma fragment DepthFrag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"

            struct AttributesD
            {
                float4 positionOS : POSITION;
            };

            struct VaryingsD
            {
                float4 positionCS : SV_POSITION;
            };

            VaryingsD DepthVert(AttributesD IN)
            {
                VaryingsD OUT;
                OUT.positionCS = TransformObjectToHClip(IN.positionOS.xyz);
                return OUT;
            }

            half4 DepthFrag(VaryingsD IN) : SV_Target { return 0; }
            ENDHLSL
        }
    }

    FallBack "Hidden/Universal Render Pipeline/FallbackError"
}
