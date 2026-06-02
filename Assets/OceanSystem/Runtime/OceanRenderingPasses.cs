using UnityEngine;
using UnityEngine.Experimental.Rendering;
using UnityEngine.Rendering;
using UnityEngine.Rendering.RenderGraphModule;
using UnityEngine.Rendering.Universal;

namespace OceanSystem
{
    public class OceanGeometryPass : ScriptableRenderPass
    {
        private readonly static ShaderTagId OceanShaderTagId = new ShaderTagId("OceanMain");
        private readonly static ProfilingSampler _profilingSampler = new ProfilingSampler("Ocean Geometry");
        private readonly OceanRendererFeature.OceanRenderingSettings _settings;
        private FilteringSettings _filteringSettings;

        public OceanGeometryPass(OceanRendererFeature.OceanRenderingSettings settings)
        {
            _settings = settings;
            _filteringSettings = new FilteringSettings(RenderQueueRange.all);
        }

        // ============================================================
        // Legacy (URP <17 / Compatibility Mode) path
        // ============================================================
        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
#if UNITY_EDITOR
            if (!OceanRendererFeature.IsRendering) return;
#endif

            CameraData cameraData = renderingData.cameraData;
            if (!OceanRendererFeature.IsCorrectCameraType(cameraData.cameraType)) return;
            Camera camera = cameraData.camera;

            DrawingSettings drawingSettings = new DrawingSettings(OceanShaderTagId,
                new SortingSettings(camera));

            CommandBuffer cmd = CommandBufferPool.Get();
            SetupGlobalKeywords(cmd);
            SetupCameraGlobals(cmd, camera);
            context.ExecuteCommandBuffer(cmd);

            using (new ProfilingScope(cmd, _profilingSampler))
            {
                context.ExecuteCommandBuffer(cmd);
                cmd.Clear();

                drawingSettings.perObjectData = PerObjectData.LightProbe;
                context.DrawRenderers(renderingData.cullResults, ref drawingSettings, ref _filteringSettings);
            }
            context.ExecuteCommandBuffer(cmd);


            CommandBufferPool.Release(cmd);
        }

        private void SetupCameraGlobals(CommandBuffer cmd, Camera cam)
        {
            cmd.SetGlobalMatrix(GlobalShaderVariables.Misc.InverseViewMatrix, cam.cameraToWorldMatrix);
            cmd.SetGlobalMatrix(GlobalShaderVariables.Misc.InverseProjectionMatrix,
                GL.GetGPUProjectionMatrix(cam.projectionMatrix, false).inverse);
        }

        private void SetupGlobalKeywords(CommandBuffer cmd)
        {
            SetGlobalKeyword(cmd, "OCEAN_TRANSPARENCY_ENABLED", _settings.transparency);
            SetGlobalKeyword(cmd, "OCEAN_UNDERWATER_ENABLED", _settings.underwaterEffect);
        }

        private void SetGlobalKeyword(CommandBuffer cmd, string keyword, bool b)
        {
            if (b)
                cmd.EnableShaderKeyword(keyword);
            else
                cmd.DisableShaderKeyword(keyword);
        }

        // ============================================================
        // RenderGraph (URP 17+) path
        // ============================================================
        private class GeometryPassData
        {
            public RendererListHandle rendererList;
            public Matrix4x4 inverseView;
            public Matrix4x4 inverseProjection;
            public bool transparency;
            public bool underwater;
        }

        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
#if UNITY_EDITOR
            if (!OceanRendererFeature.IsRendering) return;
#endif
            UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();
            if (!OceanRendererFeature.IsCorrectCameraType(cameraData.cameraType)) return;
            UniversalResourceData resourceData = frameData.Get<UniversalResourceData>();
            UniversalRenderingData renderingData = frameData.Get<UniversalRenderingData>();

            using (var builder = renderGraph.AddRasterRenderPass<GeometryPassData>(
                "Ocean Geometry", out var passData, _profilingSampler))
            {
                Camera cam = cameraData.camera;
                passData.inverseView = cam.cameraToWorldMatrix;
                passData.inverseProjection = GL.GetGPUProjectionMatrix(cam.projectionMatrix, false).inverse;
                passData.transparency = _settings.transparency;
                passData.underwater = _settings.underwaterEffect;

                DrawingSettings drawingSettings = new DrawingSettings(OceanShaderTagId, new SortingSettings(cam))
                {
                    perObjectData = PerObjectData.LightProbe
                };
                RendererListParams rlParams = new RendererListParams(
                    renderingData.cullResults, drawingSettings, _filteringSettings);
                passData.rendererList = renderGraph.CreateRendererList(rlParams);
                builder.UseRendererList(passData.rendererList);

                builder.SetRenderAttachment(resourceData.activeColorTexture, 0, AccessFlags.Write);
                builder.SetRenderAttachmentDepth(resourceData.activeDepthTexture, AccessFlags.ReadWrite);

                // Globals must propagate to subsequent passes.
                builder.AllowGlobalStateModification(true);
                // Always draw the ocean — never let RG cull this pass.
                builder.AllowPassCulling(false);

                builder.SetRenderFunc(static (GeometryPassData data, RasterGraphContext context) =>
                {
                    RasterCommandBuffer cmd = context.cmd;
                    if (data.transparency) cmd.EnableShaderKeyword("OCEAN_TRANSPARENCY_ENABLED");
                    else cmd.DisableShaderKeyword("OCEAN_TRANSPARENCY_ENABLED");
                    if (data.underwater) cmd.EnableShaderKeyword("OCEAN_UNDERWATER_ENABLED");
                    else cmd.DisableShaderKeyword("OCEAN_UNDERWATER_ENABLED");
                    cmd.SetGlobalMatrix(GlobalShaderVariables.Misc.InverseViewMatrix, data.inverseView);
                    cmd.SetGlobalMatrix(GlobalShaderVariables.Misc.InverseProjectionMatrix, data.inverseProjection);
                    cmd.DrawRendererList(data.rendererList);
                });
            }
        }
    }

    public class OceanUnderwaterEffectPass : ScriptableRenderPass
    {
        private readonly OceanRendererFeature.OceanRenderingSettings _settings;
        private readonly Material _underwaterEffectMaterial;
        private RenderTargetIdentifier _submergenceTarget;
        // URP 14+: ConfigureTarget only takes RTHandle. The temp RT still lives behind
        // a property ID (cmd.GetTemporaryRT / cmd.ReleaseTemporaryRT), but we need a
        // stable RTHandle wrapper for ConfigureTarget. RTHandles.Alloc(RenderTargetIdentifier)
        // creates a metadata-only wrapper — no GPU allocation — and the identifier itself
        // is invariant (the property ID never changes), so a single wrapper is reusable
        // forever. We allocate it lazily on first Configure and release it on disposal.
        private RTHandle _submergenceTargetRTH;
        private static readonly int _submergenceTargetID = Shader.PropertyToID("SubmergenceTarget");


        public OceanUnderwaterEffectPass(OceanRendererFeature.OceanRenderingSettings settings)
        {
            _settings = settings;
            _underwaterEffectMaterial = new Material(Shader.Find("Hidden/Ocean/UnderwaterEffect"));
        }

        // ============================================================
        // Legacy (URP <17 / Compatibility Mode) path
        // ============================================================
        public override void Configure(CommandBuffer cmd, RenderTextureDescriptor cameraTextureDescriptor)
        {
            cmd.GetTemporaryRT(_submergenceTargetID, 32, 32, 0, FilterMode.Bilinear, RenderTextureFormat.R8, RenderTextureReadWrite.Linear, 1);
            _submergenceTarget = new RenderTargetIdentifier(_submergenceTargetID);
            if (_submergenceTargetRTH == null)
                _submergenceTargetRTH = RTHandles.Alloc(_submergenceTarget, name: "SubmergenceTarget");
            ConfigureTarget(_submergenceTargetRTH);
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
#if UNITY_EDITOR
            if (!OceanRendererFeature.IsRendering) return;
#endif

            CameraData cameraData = renderingData.cameraData;
            if (!OceanRendererFeature.IsCorrectCameraType(cameraData.cameraType)) return;

            if (_settings.underwaterEffect)
            {
                CommandBuffer cmd = CommandBufferPool.Get("Underwater Effect");
                RenderingUtils.DrawPreceduralFullscreenQuad(cmd, _submergenceTarget,
                    RenderBufferLoadAction.DontCare, _underwaterEffectMaterial, 0);
                cmd.SetGlobalTexture(GlobalShaderVariables.Misc.SubmergenceTexture, _submergenceTargetID);

                // URP 13+: cameraColorTarget (RenderTargetIdentifier) was deprecated in
                // favour of cameraColorTargetHandle (RTHandle). RTHandle has an implicit
                // operator to RenderTargetIdentifier so RenderingUtils' signature still fits.
                RenderingUtils.DrawPreceduralFullscreenQuad(cmd, cameraData.renderer.cameraColorTargetHandle,
                    RenderBufferLoadAction.Load, _underwaterEffectMaterial, 1);
                context.ExecuteCommandBuffer(cmd);
                CommandBufferPool.Release(cmd);
            }
        }

        // URP 12+: ScriptableRenderPass.FrameCleanup was renamed OnCameraCleanup.
        // Same signature, same semantics — called once per camera after the pass
        // finishes, where we release the temp RT we allocated in Configure.
        public override void OnCameraCleanup(CommandBuffer cmd)
        {
            cmd.ReleaseTemporaryRT(_submergenceTargetID);
        }

        // ============================================================
        // RenderGraph (URP 17+) path
        // ============================================================
        // Implemented as a single unsafe pass so we can keep the legacy
        // RenderingUtils.DrawPreceduralFullscreenQuad call — the underwater
        // material's vertex shader expects that procedural-quad layout and
        // mis-renders (diagonal band artefact) under Blitter.BlitTexture.
        private class UnderwaterPassData
        {
            public Material material;
            public TextureHandle submergence;
            public TextureHandle cameraColor;
            public int submergenceID;
            public int submergenceGlobalID;
        }

        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
#if UNITY_EDITOR
            if (!OceanRendererFeature.IsRendering) return;
#endif
            if (!_settings.underwaterEffect) return;

            UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();
            if (!OceanRendererFeature.IsCorrectCameraType(cameraData.cameraType)) return;
            UniversalResourceData resourceData = frameData.Get<UniversalResourceData>();

            // Transient submergence target. RG owns the GPU allocation; we just
            // pass the handle through and read it as a RenderTargetIdentifier in
            // the render func via TextureHandle's implicit operator.
            TextureDesc submergenceDesc = new TextureDesc(32, 32)
            {
                colorFormat = GraphicsFormat.R8_UNorm,
                name = "SubmergenceTarget",
                filterMode = FilterMode.Bilinear,
                wrapMode = TextureWrapMode.Clamp,
                clearBuffer = false
            };
            TextureHandle submergence = renderGraph.CreateTexture(submergenceDesc);
            TextureHandle cameraColor = resourceData.activeColorTexture;

            using (var builder = renderGraph.AddUnsafePass<UnderwaterPassData>(
                "Ocean Underwater Effect", out var passData))
            {
                passData.material = _underwaterEffectMaterial;
                passData.submergence = submergence;
                passData.cameraColor = cameraColor;
                passData.submergenceID = _submergenceTargetID;
                passData.submergenceGlobalID = GlobalShaderVariables.Misc.SubmergenceTexture;

                builder.UseTexture(submergence, AccessFlags.ReadWrite);
                builder.UseTexture(cameraColor, AccessFlags.ReadWrite);
                builder.AllowGlobalStateModification(true);
                builder.AllowPassCulling(false);

                builder.SetRenderFunc(static (UnderwaterPassData data, UnsafeGraphContext context) =>
                {
                    CommandBuffer cmd = CommandBufferHelpers.GetNativeCommandBuffer(context.cmd);
                    // Submergence target — pass 0 of the underwater material.
                    RenderingUtils.DrawPreceduralFullscreenQuad(cmd, data.submergence,
                        RenderBufferLoadAction.DontCare, data.material, 0);
                    cmd.SetGlobalTexture(data.submergenceGlobalID, data.submergence);
                    // Camera color blend — pass 1 of the underwater material.
                    RenderingUtils.DrawPreceduralFullscreenQuad(cmd, data.cameraColor,
                        RenderBufferLoadAction.Load, data.material, 1);
                });
            }
        }
    }

    public class OceanSkyMapPass : ScriptableRenderPass
    {
        private OceanRendererFeature.OceanRenderingSettings _settings;


        private readonly Material _skyMapMaterial;
        private RenderTexture _skyMap;
        // URP 14+: ConfigureTarget only takes RTHandle. _skyMap is a concrete RenderTexture
        // we own and recreate when the resolution changes, so the RTHandle wrapper has to be
        // re-allocated each time we recreate the underlying RT. Released alongside it in
        // CreateSkyMapTexture.
        private RTHandle _skyMapRTH;
        private bool _skyMapRendered;
        private bool NeedToRenderSkyMap => _settings.updateSkyMap || !_skyMapRendered;

        public OceanSkyMapPass(OceanRendererFeature.OceanRenderingSettings settings)
        {
            this._settings = settings;
            _skyMapMaterial = new Material(Shader.Find("Hidden/Ocean/StereographicSky"));
        }

        // ============================================================
        // Legacy (URP <17 / Compatibility Mode) path
        // ============================================================
        public override void Configure(CommandBuffer cmd, RenderTextureDescriptor cameraTextureDescriptor)
        {
            if (NeedToRenderSkyMap)
            {
                CreateSkyMapTexture();
                ConfigureTarget(_skyMapRTH);
            }
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
#if UNITY_EDITOR
            if (!OceanRendererFeature.IsRendering) return;
#endif

            if (NeedToRenderSkyMap)
            {
                CameraData cameraData = renderingData.cameraData;
                if (!OceanRendererFeature.IsCorrectCameraType(cameraData.cameraType)) return;

                CommandBuffer cmd = CommandBufferPool.Get("Ocean Sky Map");
                RenderSkyMap(cmd);
                context.ExecuteCommandBuffer(cmd);
                CommandBufferPool.Release(cmd);
            }
        }

        private void RenderSkyMap(CommandBuffer cmd)
        {
            // URP 14+: ScriptableRenderPass.Blit was rewritten to only accept RTHandle
            // parameters. _skyMap is a raw RenderTexture (we own the allocation), so
            // we bypass the base helper and call CommandBuffer.Blit directly — that
            // overload still takes (Texture src, RenderTexture dst, Material, pass).
            cmd.Blit((Texture)null, _skyMap, _skyMapMaterial, 0);
            cmd.SetGlobalTexture(GlobalShaderVariables.Misc.SkyMap, _skyMap);
            _skyMapRendered = true;
        }

        private void CreateSkyMapTexture()
        {
            if (_skyMap == null || _skyMap.height != _settings.skyMapResolution)
            {
                if (_skyMap != null)
                    _skyMap.Release();
                // Release the prior RTHandle wrapper too — the underlying RT is being
                // replaced, so the old wrapper points at a stale handle.
                if (_skyMapRTH != null)
                {
                    RTHandles.Release(_skyMapRTH);
                    _skyMapRTH = null;
                }
                _skyMap = new RenderTexture(_settings.skyMapResolution, _settings.skyMapResolution, 0,
                    RenderTextureFormat.DefaultHDR, RenderTextureReadWrite.Linear)
                {
                    name = "SkyMap",
                    antiAliasing = 1,
                    wrapMode = TextureWrapMode.Clamp,
                    useMipMap = true,
                    autoGenerateMips = true,
                    filterMode = FilterMode.Trilinear,
                    anisoLevel = 9
                };
                _skyMap.Create();
            }
            if (_skyMapRTH == null)
                _skyMapRTH = RTHandles.Alloc(_skyMap);
        }

        // ============================================================
        // RenderGraph (URP 17+) path
        // ============================================================
        // Implemented as an unsafe pass because:
        //  (1) The legacy stereographic-sky shader expects CommandBuffer.Blit's
        //      vertex layout, not Blitter.BlitTexture's fullscreen triangle.
        //  (2) A raster pass cannot both SetRenderAttachment(_skyMap) AND
        //      SetGlobalTexture(_skyMap) in one pass (RG forbids same-frame
        //      read+write of the same texture as attachment + bound texture).
        //  Using cmd.Blit + cmd.SetGlobalTexture inside an unsafe pass sidesteps
        //  both restrictions.
        private class SkyMapPassData
        {
            public Material material;
            public RenderTexture target;
            public int globalTexID;
        }

        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
#if UNITY_EDITOR
            if (!OceanRendererFeature.IsRendering) return;
#endif
            if (!NeedToRenderSkyMap) return;

            UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();
            if (!OceanRendererFeature.IsCorrectCameraType(cameraData.cameraType)) return;

            CreateSkyMapTexture();

            using (var builder = renderGraph.AddUnsafePass<SkyMapPassData>(
                "Ocean Sky Map", out var passData))
            {
                passData.material = _skyMapMaterial;
                passData.target = _skyMap;
                passData.globalTexID = GlobalShaderVariables.Misc.SkyMap;

                builder.AllowGlobalStateModification(true);
                builder.AllowPassCulling(false);

                builder.SetRenderFunc(static (SkyMapPassData data, UnsafeGraphContext context) =>
                {
                    CommandBuffer cmd = CommandBufferHelpers.GetNativeCommandBuffer(context.cmd);
                    cmd.Blit((Texture)null, data.target, data.material, 0);
                    cmd.SetGlobalTexture(data.globalTexID, data.target);
                });
            }
            _skyMapRendered = true;
        }
    }
}
