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

            // ---- RTS hover uniforms ----
            //
            // Set by TileSelector (Shader.SetGlobalXxx) on cursor move:
            //   _HoverPrimalIdx ∈ [0, chunkCellCount); -1 = nothing hovered
            //   _HoverChunkQR    = the hovered chunk's (q, r) axial coord (.zw unused)
            // Set by GameModeController on Tab toggle:
            //   _HoverColor      = (rgb, a) tint blended over the hovered cell
            // Set per renderer (MaterialPropertyBlock from ChunkSurface):
            //   _ChunkQR         = THIS chunk's (q, r) coord; only chunks whose coord
            //                      matches _HoverChunkQR participate in the hover tint
            //                      (otherwise the same local index would tint cells in
            //                      every chunk simultaneously).
            //
            // Living at global scope (outside UnityPerMaterial CBUFFER) so MPB overrides
            // work cleanly and SetGlobal* targets a single shared backing. SRP Batcher
            // compatibility is unaffected by globals; per-renderer MPB on _ChunkQR will
            // break batching for terrain chunks, which is acceptable cost (~50 chunks
            // around the camera at typical renderRadius).
            //
            // (Outline rendering does NOT live in this shader — see ChunkSurface
            // BuildOutlineMesh for the per-chunk line-topology overlay.)
            float  _HoverPrimalIdx;
            float4 _HoverChunkQR;
            float4 _HoverColor;
            float4 _ChunkQR;

            struct Attributes
            {
                float4 positionOS : POSITION;
                float3 normalOS   : NORMAL;
                float4 color      : COLOR;       // per-vertex tile-property tint (legacy / fallback)
                // Nearest-neighbour tile shading: BuildMesh unshares vertices per triangle and
                // packs (barycentric, cornerA_rgb, cornerB_rgb, cornerC_rgb) onto every vertex
                // of the triangle. All three vertices of a triangle output the SAME three
                // corner colours so interpolation leaves the tuple constant across the
                // triangle; the barycentric, however, varies per vertex and interpolates to
                // the actual barycentric weight at each fragment.
                float3 barycentric : TEXCOORD0;
                float3 cornerA     : TEXCOORD1;
                float3 cornerB     : TEXCOORD2;
                float3 cornerC     : TEXCOORD3;
                // Per-triangle primal cell indices of the three corners. Same pattern
                // as the corner colours — every vertex of a triangle carries the same
                // triple, so interpolation leaves it constant across the fragment area.
                // Used for the shader-side hover tint and per-cell outline detection
                // (a fragment is "on a cell boundary" iff its two largest barycentric
                // weights belong to corners with different primal indices). -1 means
                // foreign corner (cross-chunk seam) — excluded from both hover and
                // outline so cross-chunk seams don't draw spurious contours.
                float3 cornerPrimal : TEXCOORD4;
            };

            struct Varyings
            {
                float4 positionCS   : SV_POSITION;
                float3 positionWS   : TEXCOORD0;
                float3 normalWS     : TEXCOORD1;
                float  fogFactor    : TEXCOORD2;
                float3 barycentric  : TEXCOORD3;
                float3 cornerA      : TEXCOORD4;
                float3 cornerB      : TEXCOORD5;
                float3 cornerC      : TEXCOORD6;
                float3 cornerPrimal : TEXCOORD7;
                float4 color        : COLOR;      // legacy pass-through (still consumed for flat tiles)
            };

            Varyings Vert(Attributes IN)
            {
                Varyings OUT;
                VertexPositionInputs pos = GetVertexPositionInputs(IN.positionOS.xyz);
                VertexNormalInputs   nrm = GetVertexNormalInputs(IN.normalOS);
                OUT.positionCS   = pos.positionCS;
                OUT.positionWS   = pos.positionWS;
                OUT.normalWS     = nrm.normalWS;
                OUT.fogFactor    = ComputeFogFactor(pos.positionCS.z);
                OUT.color        = IN.color;
                OUT.barycentric  = IN.barycentric;
                OUT.cornerA      = IN.cornerA;
                OUT.cornerB      = IN.cornerB;
                OUT.cornerC      = IN.cornerC;
                OUT.cornerPrimal = IN.cornerPrimal;
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

                // 2b) Tile-property tint, nearest-neighbour across each triangle. The mesh
                //     ships each triangle as 3 unshared vertices carrying the same triple of
                //     corner colours (cornerA/B/C) plus barycentric identifiers that
                //     interpolate to per-pixel weights. We pick the colour whose weight is
                //     largest — i.e. the Voronoi region of the nearest primal corner —
                //     giving sharp uniformly-coloured primal-cell regions instead of the
                //     soft gradient that linear interpolation produced before. Multiplies
                //     AFTER biome / slope so a painted road still reads slightly grassy on
                //     grass and slightly rocky on a cliff.
                //
                //     Also captures which corner index won (0/1/2) so the hover tint
                //     below uses the same nearest-neighbour decision.
                float3 b = IN.barycentric;
                int nearest;
                float3 tileTint;
                if (b.x >= b.y && b.x >= b.z)      { tileTint = IN.cornerA; nearest = 0; }
                else if (b.y >= b.z)               { tileTint = IN.cornerB; nearest = 1; }
                else                               { tileTint = IN.cornerC; nearest = 2; }
                col *= tileTint;

                // 2c) Hover tint — pixel-exact across the rendered cell footprint.
                //     Fires only when (a) this fragment's nearest corner belongs to
                //     the hovered primal cell index AND (b) this chunk's coord matches
                //     the hovered chunk. The second test exists because primal indices
                //     are local to each chunk: the same index 12 can mean different
                //     cells in different chunks, so without the chunk gate we'd tint
                //     cell 12 in every chunk on screen at once.
                //
                //     Outlines do NOT live in this shader anymore — barycentric Voronoi
                //     bisectors form Y-shapes meeting at triangle centroids, not the
                //     wedge perimeters players perceive as tile boundaries. The per-
                //     chunk outline overlay built by ChunkSurface (line topology mesh
                //     of V → M segments) handles outlines instead; it coincides
                //     exactly with the wedge boundaries that BuildMesh emits.
                float idxN = (nearest == 0) ? IN.cornerPrimal.x
                          : (nearest == 1) ? IN.cornerPrimal.y
                                           : IN.cornerPrimal.z;
                bool chunkMatch = abs(_ChunkQR.x - _HoverChunkQR.x) < 0.5
                               && abs(_ChunkQR.y - _HoverChunkQR.y) < 0.5;
                bool idxMatch   = idxN >= 0.0 && abs(idxN - _HoverPrimalIdx) < 0.5;
                bool hoverHere  = _HoverPrimalIdx >= 0.0 && chunkMatch && idxMatch;
                col = lerp(col, _HoverColor.rgb, hoverHere ? _HoverColor.a : 0.0);

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
