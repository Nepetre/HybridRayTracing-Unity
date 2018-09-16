using System.Collections;
using System.Collections.Generic;
using UnityEngine;


[ExecuteInEditMode]
[RequireComponent(typeof(Camera))]
public class HybridRayTracer : MonoBehaviour
{

    #region Private variables
    private Camera camA = null;
    private GameObject camAGO = null;

    private uint[] initHead;

    private int Width { get { return GetComponent<Camera>().pixelWidth; } }
    private int Height { get { return GetComponent<Camera>().pixelHeight; } }

    private RenderTexture prevFrame = null;

    private Vector3 oldCamPos = Vector3.zero;
    private Vector3 oldCamDir = Vector3.zero;
    private int framecountOld = 0;

    private float pixelStrideZCutoff = 0.0f;
    private float pixelZSizeOffset = 0.0f;
    private float maxRayDistance = 1.0f;

    private int maxIntIteration = 0;

    private Shader _backfaceDepthShader;
    private Material _ssrMaterial;
    private Material _blurMaterial;
    private Material _combinerMaterial;
    private RenderTexture _backFaceDepthTexture;
    private Camera _camera = null;
    private Camera _backFaceCamera;

    private bool hasMoved = false;
    #endregion

    #region Public variables
    [Header("Downsample:")]
    [Range(0.0f, 8.0f)]
    public int backfaceDepthDownsample = 0;
    [Range(0.0f, 8.0f)]
    public int ssrDownsample = 0;
    [Range(1.0f, 36.0f)]
    public int perPixelLinkedListDepth = 1;

    [Header("Raycast:")]
    public bool Refinement = true;
    [Range(0.0f, 50.0f)]
    public int rayPP = 1;
    [Range(0.0f, 20.0f)]
    public int numBounces = 4;
    [Range(1.0f, 2200.0f)]
    public int iterations = 400;
    public bool maxIteration = false;
    public int maxIterationDebug = 0;
    [Range(1.0f, 100.0f)]
    public int pixelStride = 1;
    #endregion
    
    private void Awake()
    {
        _camera = GetComponent<Camera>();
        _camera.renderingPath = RenderingPath.DeferredShading;
        //_camera.depthTextureMode |= DepthTextureMode.Depth;
        if (camAGO != null)
        {
            DestroyImmediate(camAGO);
        }
        camAGO = new GameObject("Per-Pixel-Linked-List-Cam")
        {
            hideFlags = HideFlags.DontSave
        };
        camAGO.transform.parent = transform;
        camAGO.transform.localPosition = Vector3.zero;
        camA = camAGO.AddComponent<Camera>();
        camA.CopyFrom(_camera);
        camA.clearFlags = CameraClearFlags.SolidColor;
        camA.enabled = false;
        //camA.depthTextureMode |= DepthTextureMode.Depth;

        initHead = new uint[Height * Width];
        for (int i = 0; i < initHead.Length; i++)
        {
            initHead[i] = 0xffffffff;
        }

        prevFrame = new RenderTexture(Width, Height, 0, RenderTextureFormat.Default);
        prevFrame.Create();

        maxIntIteration = (int)(Mathf.Sqrt(Width * Width + Height * Height));
    }

    private void OnDestroy()
    {
        DestroyImmediate(camAGO);
    }

    private void OnPreRender()
    {
        //cam.cullingMask = 0;
    }

    void OnPreCull()
    {
        
        int downsampleBackFaceDepth = backfaceDepthDownsample + 1;
        int width = _camera.pixelWidth / downsampleBackFaceDepth;
        int height = _camera.pixelHeight / downsampleBackFaceDepth;
        _backFaceDepthTexture = RenderTexture.GetTemporary(width,
                                                           height,
                                                           16,
                                                           RenderTextureFormat.RFloat);

        if (_backFaceCamera == null)
        {
            GameObject cameraGameObject = new GameObject("BackFaceDepthCamera");
            cameraGameObject.hideFlags = HideFlags.DontSave;
            cameraGameObject.transform.parent = transform;
            cameraGameObject.transform.localPosition = Vector3.zero;
            _backFaceCamera = cameraGameObject.AddComponent<Camera>();
        }

        _backFaceCamera.CopyFrom(_camera);
        _backFaceCamera.renderingPath = RenderingPath.Forward;
        _backFaceCamera.enabled = false;
        _backFaceCamera.SetReplacementShader(Shader.Find("Hidden/BackfaceDepth"), null);
        _backFaceCamera.backgroundColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
        _backFaceCamera.clearFlags = CameraClearFlags.SolidColor;
        //_backFaceCamera.cullingMask = LayerMask.GetMask( "Everything");

        _backFaceCamera.targetTexture = _backFaceDepthTexture;
        _backFaceCamera.Render();
        
    }

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        #region Per-Pixel Linked list
        

        Vector3 camPos = GetComponent<Camera>().transform.position;
        Vector3 camDir = GetComponent<Camera>().transform.forward;
        Vector3 camUp = GetComponent<Camera>().transform.up;
        Vector3 camRight = GetComponent<Camera>().transform.right;

        int it = maxIteration ? maxIntIteration : iterations;

        maxIterationDebug = it;

        int framecount = Time.frameCount;
        int realFramecount = framecount;

      
        if (camPos != oldCamPos || camDir != oldCamDir || hasMoved)
        {
            framecountOld = framecount;
            oldCamPos = camPos;
            oldCamDir = camDir;
        }
        
        hasMoved = false;
        framecount -= framecountOld;
        framecount += rayPP;


        ComputeBuffer nodeBuffer = new ComputeBuffer((Width * Height) * perPixelLinkedListDepth, 20 * sizeof(float) + sizeof(uint), ComputeBufferType.Default);
        nodeBuffer.SetCounterValue(0);

        ComputeBuffer head = new ComputeBuffer((Width * Height), sizeof(uint), ComputeBufferType.Raw);
        head.SetData(initHead);


        ComputeBuffer counter = new ComputeBuffer(1, sizeof(uint), ComputeBufferType.Counter);
        counter.SetCounterValue(0);

        Graphics.ClearRandomWriteTargets();

        Graphics.SetRandomWriteTarget(1, nodeBuffer);
        Graphics.SetRandomWriteTarget(2, head);
        Graphics.SetRandomWriteTarget(3, counter);

        Shader.SetGlobalBuffer("list", nodeBuffer);
        Shader.SetGlobalBuffer("head", head);



        Shader.SetGlobalInt("_RPP", rayPP);
        Shader.SetGlobalInt("width", Width);
        Shader.SetGlobalInt("_RayBounce", numBounces);
        
        Shader.SetGlobalTexture("_BackFaceDepthTex", _backFaceDepthTexture);
        Shader.SetGlobalVector("_OneDividedByRenderBufferSize", new Vector4(1.0f / Width, 1.0f / Height, 0.0f, 0.0f));
        Shader.SetGlobalFloat("_PixelZSize", pixelZSizeOffset);
        Shader.SetGlobalFloat("_MaxRayDistance", maxRayDistance);
        
        Shader.SetGlobalFloat("_PixelStride", pixelStride);
        Shader.SetGlobalFloat("_PixelStrideZCuttoff", pixelStrideZCutoff);
        Shader.SetGlobalFloat("_PixelZSize", pixelZSizeOffset);

        Shader.SetGlobalFloat("_Iterations", it);
        Shader.SetGlobalFloat("_MaxRayDistance", maxRayDistance);

        int downsampleSSR = ssrDownsample + 1;
        int width = _camera.pixelWidth / downsampleSSR;
        int height = _camera.pixelHeight / downsampleSSR;
        Matrix4x4 trs = Matrix4x4.TRS(new Vector3(0.5f, 0.5f, 0.0f), Quaternion.identity, new Vector3(0.5f, 0.5f, 1.0f));
        Matrix4x4 scrScale = Matrix4x4.Scale(new Vector3(width, height, 1.0f));
        Matrix4x4 projection = _camera.projectionMatrix;

        Matrix4x4 m = scrScale * trs * projection;

        Shader.SetGlobalVector("_RenderBufferSize", new Vector4(width, height, 0.0f, 0.0f));
        Shader.SetGlobalVector("_OneDividedByRenderBufferSize", new Vector4(1.0f / width, 1.0f / height, 0.0f, 0.0f));
        Shader.SetGlobalMatrix("_CameraProjectionMatrix", m);
        Shader.SetGlobalMatrix("_CameraInverseProjectionMatrix", projection.inverse);
        Shader.SetGlobalMatrix("_NormalMatrix", _camera.worldToCameraMatrix);
        
        Shader.SetGlobalVector("resolution", new Vector2(Width, Height));
        
        Shader.SetGlobalInt("_Refinement", Refinement == true ? 1 : 0);
        Shader.SetGlobalTexture("_PrevFrame", prevFrame);
        Shader.SetGlobalInt("_FrameCount", framecount);
        Shader.SetGlobalInt("_RealFrameCount", realFramecount);

        Shader.SetGlobalVector("camPos", new Vector4(camPos.x, camPos.y, camPos.z, 1));
        Shader.SetGlobalVector("camDir", new Vector4(camDir.x, camDir.y, camDir.z, 1));
        Shader.SetGlobalVector("camUp", new Vector4(camUp.x, camUp.y, camUp.z, 1));
        Shader.SetGlobalVector("camRight", new Vector4(camRight.x, camRight.y, camRight.z, 1));

        
        float nearClip = camA.nearClipPlane;
        Matrix4x4 view = camA.worldToCameraMatrix;
        Matrix4x4 p = GL.GetGPUProjectionMatrix(camA.projectionMatrix, false);
        Matrix4x4 VP = p * view;


        Shader.SetGlobalFloat("_nearClip", nearClip);
        Shader.SetGlobalMatrix("_V", view);
        Shader.SetGlobalMatrix("_P", p);
        Shader.SetGlobalMatrix("_VP", VP);
        Shader.SetGlobalMatrix("_InvVP", VP.inverse);
        Shader.SetGlobalMatrix("_InvP", p.inverse);
        Shader.SetGlobalMatrix("_InvV", view.inverse);
        

        if (camA == null)
        {
            GameObject newCamA = new GameObject("New Per-Pixel-Linked-List-Cam");
            newCamA.transform.parent = transform;
            newCamA.transform.localPosition = Vector3.zero;
            newCamA.hideFlags = HideFlags.DontSave;
            camA = newCamA.AddComponent<Camera>();
        }

        camA.CopyFrom(_camera);
        camA.renderingPath = RenderingPath.Forward;
        camA.enabled = false;
        camA.SetReplacementShader(Shader.Find("Hidden/PerPixelLinkedList"), null);
        camA.backgroundColor = new Color(1.0f, 1.0f, 1.0f, 1.0f);
        camA.clearFlags = CameraClearFlags.SolidColor;
        
        camA.Render();
        #endregion



        #region Raytracing

        Graphics.ClearRandomWriteTargets();

        RenderTexture testTex = RenderTexture.GetTemporary(Width, Height, 0, RenderTextureFormat.Default);

        Material mat = new Material(Shader.Find("Hidden/HybridRayTracing"));
        mat.SetTexture("_PrevFrame", prevFrame);
        //mat.hideFlags = HideFlags.HideAndDontSave;
        Graphics.Blit(source, testTex, mat);

        Graphics.Blit(testTex, prevFrame);

        Graphics.Blit(testTex, destination);

        #endregion

        #region Release Buffers
        nodeBuffer.Release();
        head.Release();
        counter.Release();
        RenderTexture.ReleaseTemporary(testTex);
        RenderTexture.ReleaseTemporary(_backFaceDepthTexture);
        #endregion



    }

    public void HasMoved()
    {
        hasMoved = true;
    }
    
}
