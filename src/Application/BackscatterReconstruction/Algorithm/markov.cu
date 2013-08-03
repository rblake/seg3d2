
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

#define MAX_GPUS 2


#include <Application/BackscatterReconstruction/Algorithm/markov.h>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>

#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_functions.h>
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_cuda.h>
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/cutil_math.h>

//#define CUDA_ENABLE_UPDATE_FORWARD_PROJECTION


#ifdef WIN32
#include <windows.h>
class Timer {
  public:
  Timer(const char *name) :
    mName(name) {
    mTotalTime.QuadPart = 0;
    mStartTime.QuadPart = 0;
  }

  ~Timer() {
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    std::cerr<<mName<<": "<<(double)mTotalTime.QuadPart / frequency.QuadPart<<std::endl;
  }


  void Start() {
    QueryPerformanceCounter(&mStartTime);
  }

  void Stop() {
    LARGE_INTEGER stopTime;
    QueryPerformanceCounter(&stopTime);
    mTotalTime.QuadPart += stopTime.QuadPart - mStartTime.QuadPart;
  }


  const char *mName;
  LARGE_INTEGER mTotalTime;
  LARGE_INTEGER mStartTime;
};

#if 0
#define TIMER_START(name)   static Timer timer(name);  timer.Start();
#define TIMER_STOP                              \
  for (int g=0; g<mGpuIds.size(); g++) {        \
    checkCudaErrors(cudaSetDevice(mGpuIds[g])); \
    checkCudaErrors(cudaDeviceSynchronize());   \
  }                                             \
  timer.Stop();
#else
#define TIMER_START(name)
#define TIMER_STOP
#endif

#else

#define TIMER_START(name)
#define TIMER_STOP

#endif



// all the memory / texture references needed for each gpu.  Texture references 
// are magical and automatically get duplicated for each gpu.

// material volumes in the collection
float4 *cudaVolumeLinearCollection[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray00[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray01[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray02[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray10[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray11[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray12[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray20[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray21[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaVolumeArray22[MAX_GPUS] = { NULL, NULL };
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures00;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures01;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures02;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures10;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures11;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures12;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures20;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures21;
texture<float4, cudaTextureType3D, cudaReadModeElementType> cudaVolumeTextures22;


// source attenuation volumes in the collection
cudaArray *cudaSourceArray00[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray01[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray02[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray10[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray11[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray12[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray20[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray21[MAX_GPUS] = { NULL, NULL };
cudaArray *cudaSourceArray22[MAX_GPUS] = { NULL, NULL };
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures00;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures01;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures02;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures10;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures11;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures12;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures20;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures21;
texture<float, cudaTextureType3D, cudaReadModeElementType> cudaSourceTextures22;


// material scatter factors
texture<float, cudaTextureType1D, cudaReadModeElementType> cudaMaterialTextures0;
texture<float, cudaTextureType1D, cudaReadModeElementType> cudaMaterialTextures1;
texture<float, cudaTextureType1D, cudaReadModeElementType> cudaMaterialTextures2;
cudaArray* cudaMaterialArray[MAX_GPUS][NUM_MATERIALS];


// source ray info
float3 *cudaSourceRayVolOrigin[MAX_GPUS] = { NULL, NULL };
float3 *cudaSourceRayVolDir[MAX_GPUS] = { NULL, NULL };
float *cudaSourceRayTMin[MAX_GPUS] = { NULL, NULL };
float *cudaSourceRayTMax[MAX_GPUS] = { NULL, NULL };
float *cudaSourceScale[MAX_GPUS] = { NULL, NULL }; // accounts for source 1/r^2 attenuation
float *cudaSourceAttenuationCollection[MAX_GPUS] = { NULL, NULL }; // output
texture<float, cudaTextureType2D, cudaReadModeElementType> cudaSourceFalloffTexture;
cudaArray *cudaFalloffArray[MAX_GPUS] = { NULL, NULL };
float *cudaProjectionAngles[MAX_GPUS] = { NULL, NULL };

// detector ray info
float3 *cudaDetectorRayWorldOrigin[MAX_GPUS] = { NULL, NULL };
float3 *cudaDetectorRayWorldDir[MAX_GPUS] = { NULL, NULL };
float3 *cudaDetectorRayVolOrigin[MAX_GPUS] = { NULL, NULL };
float3 *cudaDetectorRayVolDir[MAX_GPUS] = { NULL, NULL };
float *cudaDetectorRayTMin[MAX_GPUS] = { NULL, NULL };
float *cudaDetectorRayTMax[MAX_GPUS] = { NULL, NULL };
float *cudaForwardProjectionCollection[MAX_GPUS] = { NULL, NULL }; // output

// info for calculating projection errors
float *cudaBaselineProjection[MAX_GPUS] = { NULL, NULL };
float *cudaForwardProjectionError[MAX_GPUS] = { NULL, NULL };


// info for optimizing which rays to cast
unsigned int *cudaRayIds[MAX_GPUS] = { NULL, NULL };
float *cudaRayPriority[MAX_GPUS] = { NULL, NULL };
float *cudaCurrentForwardProjection[MAX_GPUS] = { NULL, NULL };


__device__ float sampleScatterFactor(int mat, float x) {
  // cell centered to node centered indexing with normalized coordinates
  x = x*((MATERIAL_ANGULAR_SAMPLES-1.0f)/MATERIAL_ANGULAR_SAMPLES) + (0.5f/MATERIAL_ANGULAR_SAMPLES);

  switch (mat) {
  case 0:
    return tex1D(cudaMaterialTextures0, x);
  case 1:
    return tex1D(cudaMaterialTextures1, x);
  case 2:
    return tex1D(cudaMaterialTextures2, x);
  }
  return 0;
}


inline __host__ __device__ float3 operator-(float3 a, int3 b)
{
  return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ float3 operator-(int3 a, float3 b)
{
  return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ float3 operator/(float3 a, int3 b)
{
  return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}


__device__ float sampleSourceAttenuation(int combo, float3 x, int3 volDim) {
  // normalize coordinate
  x = x / (volDim-make_int3(1));

  // cell centered to node centered indexing with normalized coordinates
  x = x*((volDim-make_float3(1.0f))/volDim) + (make_float3(0.5f)/volDim);

  switch (combo) {
  case 0:  return tex3D(cudaSourceTextures00, x.x, x.y, x.z);
  case 1:  return tex3D(cudaSourceTextures01, x.x, x.y, x.z);
  case 2:  return tex3D(cudaSourceTextures02, x.x, x.y, x.z);
  case 3:  return tex3D(cudaSourceTextures10, x.x, x.y, x.z);
  case 4:  return tex3D(cudaSourceTextures11, x.x, x.y, x.z);
  case 5:  return tex3D(cudaSourceTextures12, x.x, x.y, x.z);
  case 6:  return tex3D(cudaSourceTextures20, x.x, x.y, x.z);
  case 7:  return tex3D(cudaSourceTextures21, x.x, x.y, x.z);
  case 8:  return tex3D(cudaSourceTextures22, x.x, x.y, x.z);
  }
  return 0; // should never get here
}


__device__ void sampleMaterialsAtPoint(int combo, float3 x, int3 volDim, int3 *mats, float3 *conc) {
  // normalize coordinate
  x = x / (volDim-make_int3(1));

  // cell centered to node centered indexing with normalized coordinates
  x = x*((volDim-make_float3(1.0f))/volDim) + (make_float3(0.5f)/volDim);

  float4 sample;
  switch (combo) {
  case 0:   sample = tex3D(cudaVolumeTextures00, x.x, x.y, x.z);  break;
  case 1:   sample = tex3D(cudaVolumeTextures01, x.x, x.y, x.z);  break;
  case 2:   sample = tex3D(cudaVolumeTextures02, x.x, x.y, x.z);  break;
  case 3:   sample = tex3D(cudaVolumeTextures10, x.x, x.y, x.z);  break;
  case 4:   sample = tex3D(cudaVolumeTextures11, x.x, x.y, x.z);  break;
  case 5:   sample = tex3D(cudaVolumeTextures12, x.x, x.y, x.z);  break;
  case 6:   sample = tex3D(cudaVolumeTextures20, x.x, x.y, x.z);  break;
  case 7:   sample = tex3D(cudaVolumeTextures21, x.x, x.y, x.z);  break;
  case 8:   sample = tex3D(cudaVolumeTextures22, x.x, x.y, x.z);  break;
  }


  // sort the concentrations so highest concentration is at .x, second is at .y
  (*mats) = make_int3(0,1,2);
  (*conc) = make_float3(sample.x, sample.y, sample.z);

  if (conc->x < conc->y) {
    float ti = mats->x;
    mats->x = mats->y;
    mats->y = ti;

    float tf = conc->x;
    conc->x = conc->y;
    conc->y = tf;
  }

  if (conc->y < conc->z) {
    float ti = mats->y;
    mats->y = mats->z;
    mats->z = ti;

    float tf = conc->y;
    conc->y = conc->z;
    conc->z = tf;
  }

  if (conc->x < conc->y) {
    float ti = mats->x;
    mats->x = mats->y;
    mats->y = ti;

    float tf = conc->x;
    conc->x = conc->y;
    conc->y = tf;
  }
}


//==================================================================================================
//==================================================================================================
//==================================================================================================
__global__ void castSourceRays(const float *cudaSourceScale,
                               const float3 *cudaSourceRayVolOrigin,
                               const float3 *cudaSourceRayVolDir,
                               const float *cudaSourceRayTMin,
                               const float *cudaSourceRayTMax,
                               int3 volDim,
                               int collectionStart,
                               int collectionEnd,
                               float3 matAtten,
                               float3 matDensity,
                               float voxelStepSize,
                               float *cudaSourceAttenuationCollection) {

  unsigned int gx = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int gy = blockIdx.y*blockDim.y + threadIdx.y;
  unsigned int gz = blockIdx.z*blockDim.z + threadIdx.z;
  if (gx>=volDim.x || gy>=volDim.y || 
      gz<volDim.z*collectionStart || gz>=volDim.z*collectionEnd)
    return;
  unsigned int rayIndex = (gz%volDim.z)*volDim.x*volDim.y + gy*volDim.x + gx;
  int combo = gz/volDim.z;


  float3 volOrigin = cudaSourceRayVolOrigin[rayIndex];
  float3 volDir = cudaSourceRayVolDir[rayIndex];
  float tmin = cudaSourceRayTMin[rayIndex];
  float tmax = cudaSourceRayTMax[rayIndex];

  if (tmin>tmax) {
    // no attenuation from the volume
    cudaSourceAttenuationCollection[combo*volDim.x*volDim.y*volDim.z + rayIndex] = cudaSourceScale[rayIndex];
    return;
  }

  float3 volP0 = volOrigin + tmin*volDir;
  float3 volP1 = volOrigin + tmax*volDir;

  float volLength = length(volP1-volP0);
  int numSteps = volLength / voxelStepSize + 1;


  float lastInterfaceT = tmin;
  int3 lastMat;
  float3 lastConc;

  float3 matLens = make_float3(0,0,0);

  for (int step=0; step<numSteps; step++) {
    float samplef = (step+0.5f) / numSteps;
    float3 volSample = lerp(volP0, volP1, samplef);

    int3 thisMat;
    float3 thisConc;
    sampleMaterialsAtPoint(combo, volSample, volDim, &thisMat, &thisConc);

    // integrate between last step and this one
    if (step > 0) {
      float interfacef;

      // interpolate an interface
      if (thisMat.x != lastMat.x) {
        interfacef = (lastConc.x-lastConc.y) / ((thisConc.x-lastConc.y) - (thisConc.y-lastConc.x));
        interfacef = clamp(interfacef, 0.0f, 1.0f);
        interfacef = lerp((step-0.5f) / numSteps, samplef, interfacef);
      }

      // fixed interface
      else {
        interfacef = (float)step / numSteps;
      }

      float interfaceT = lerp(tmin, tmax, interfacef);

      // add material length before current interface
      float len = interfaceT - lastInterfaceT;
      switch (lastMat.x) {
      case 0:  matLens.x += len;  break;
      case 1:  matLens.y += len;  break;
      case 2:  matLens.z += len;  break;
      }

      lastInterfaceT = interfaceT;
    }

    lastMat = thisMat;
    lastConc = thisConc;
  }

  // add final material length
  float len = tmax - lastInterfaceT;
  switch (lastMat.x) {
  case 0:  matLens.x += len;  break;
  case 1:  matLens.y += len;  break;
  case 2:  matLens.z += len;  break;
  }


  // intergrate attenuations
  cudaSourceAttenuationCollection[combo*volDim.x*volDim.y*volDim.z + rayIndex] =
    cudaSourceScale[rayIndex] * 
    exp(-matLens.x * matDensity.x * matAtten.x +
        -matLens.y * matDensity.y * matAtten.y +
        -matLens.z * matDensity.z * matAtten.z);
}
                               


void MarkovContext::CudaComputeSourceAttenuation(int collectionSize) const {
  TIMER_START("CudaComputeSourceAttenuation()");

  int3 volDim = make_int3(mGeometry.GetVolumeNodeSamplesX(),
                          mGeometry.GetVolumeNodeSamplesY(),
                          mGeometry.GetVolumeNodeSamplesZ());

  float3 matAtten = make_float3(mMaterials[0].GetMassAttenuationCoefficient(),
                                mMaterials[1].GetMassAttenuationCoefficient(),
                                mMaterials[2].GetMassAttenuationCoefficient());
  float3 matDensity = make_float3(mMaterials[0].GetDensity(),
                                  mMaterials[1].GetDensity(),
                                  mMaterials[2].GetDensity());
                                

  dim3 dimBlock(8, 8, 1);
  dim3 dimGrid(1+(volDim.x-1) / dimBlock.x, 
               1+(volDim.y-1) / dimBlock.y, 
               1+(collectionSize*volDim.z-1) / dimBlock.z);

  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, collectionSize, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    castSourceRays<<<dimGrid, dimBlock>>>(cudaSourceScale[g],
                                          cudaSourceRayVolOrigin[g],
                                          cudaSourceRayVolDir[g],
                                          cudaSourceRayTMin[g],
                                          cudaSourceRayTMax[g],
                                          volDim,
                                          collectionStart,
                                          collectionEnd,
                                          matAtten,
                                          matDensity,
                                          mVoxelStepSize,
                                          cudaSourceAttenuationCollection[g]);

    for (int c=collectionStart; c<collectionEnd; c++) {
      cudaArray **cudaArray = NULL;
      switch (c) {
      case 0:
        cudaArray = &cudaSourceArray00[g];
        break;
      case 1:
        cudaArray = &cudaSourceArray01[g];
        break;
      case 2:
        cudaArray = &cudaSourceArray02[g];
        break;
      case 3:
        cudaArray = &cudaSourceArray10[g];
        break;
      case 4:
        cudaArray = &cudaSourceArray11[g];
        break;
      case 5:
        cudaArray = &cudaSourceArray12[g];
        break;
      case 6:
        cudaArray = &cudaSourceArray20[g];
        break;
      case 7:
        cudaArray = &cudaSourceArray21[g];
        break;
      case 8:
        cudaArray = &cudaSourceArray22[g];
        break;
      }


      cudaExtent volumeSize = make_cudaExtent(mGeometry.GetVolumeNodeSamplesX(),
                                              mGeometry.GetVolumeNodeSamplesY(),
                                              mGeometry.GetVolumeNodeSamplesZ());
      // copy data to 3D array
      cudaMemcpy3DParms copyParams = {0};
      copyParams.srcPtr   = make_cudaPitchedPtr((void *)&cudaSourceAttenuationCollection[g][mGeometry.GetTotalVolumeNodeSamples() * c],
                                                volumeSize.width*sizeof(float), volumeSize.width, volumeSize.height);
      copyParams.dstArray = *cudaArray;
      copyParams.extent   = volumeSize;
      copyParams.kind     = cudaMemcpyDeviceToDevice;
      checkCudaErrors(cudaMemcpy3DAsync(&copyParams));
    }
  }

  TIMER_STOP;
}

void MarkovContext::CudaGetSourceAttenuation(int collectionSize,
                                             vector< vector<float> > &sourceAttenuationCollection) const {
  // copy results back 
  vector<float> outputData(collectionSize * mGeometry.GetTotalVolumeNodeSamples());
  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, collectionSize, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    if ((collectionEnd-collectionStart) > 0) {
      checkCudaErrors(cudaMemcpy(&outputData[collectionStart * mGeometry.GetTotalVolumeNodeSamples()],
                                 &cudaSourceAttenuationCollection[g][collectionStart * mGeometry.GetTotalVolumeNodeSamples()],
                                 (collectionEnd-collectionStart) * mGeometry.GetTotalVolumeNodeSamples()*sizeof(float),
                                 cudaMemcpyDeviceToHost));
    }
  }

  sourceAttenuationCollection.resize(collectionSize);
  for (int c=0; c<collectionSize; c++) {
    sourceAttenuationCollection[c].resize(mCurrentVolumeSourceAttenuation.size());
    for (int i=0; i<mCurrentVolumeSourceAttenuation.size(); i++)
      sourceAttenuationCollection[c][i] = outputData[c*mGeometry.GetTotalVolumeNodeSamples() + i];
  }
}




//==================================================================================================
//==================================================================================================
//==================================================================================================

__global__ void prioritizeDetectorRays(const float3 *cudaDetectorRayWorldOrigin,
                                       const float3 *cudaDetectorRayWorldDir,
                                       const float *cudaDetectorRayTMin,
                                       const float *cudaDetectorRayTMax,
                                       int numRays,
                                       float3 coneOrigin,
                                       float3 coneDir,
                                       float coneCosTheta,
                                       float coneMinDist,
                                       unsigned int *cudaRayIds,
                                       float *cudaRayPriority) {
  unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx >= numRays)
    return;

  float3 rayOrigin = cudaDetectorRayWorldOrigin[idx];
  float3 rayDir = cudaDetectorRayWorldDir[idx];


  float AdD = dot(coneDir, rayDir);
  float cosSqr = (coneCosTheta-COSTHETA_EPS)*(coneCosTheta-COSTHETA_EPS);
  float3 E = rayOrigin - coneOrigin;
  float AdE = dot(coneDir, E);
  float DdE = dot(rayDir, E);
  float EdE = dot(E, E);
  float c2 = AdD*AdD - cosSqr;
  float c1 = AdD*AdE - cosSqr*DdE;
  float c0 = AdE*AdE - cosSqr*EdE;
  float dp;

  float3 point;
  bool hit = false;

  // Solve the quadratic.  Keep only those X for which Dot(A,X-V) >= 0.
  if (fabsf(c2) >= 1e-4) {
    // c2 != 0
    float discr = c1*c1 - c0*c2;
    if (discr > 1e-4) {
      // Q(t) = 0 has two distinct real-valued roots.  However, one or
      // both of them might intersect the portion of the double-sided
      // cone "behind" the vertex.  We are interested only in those
      // intersections "in front" of the vertex.
      float root = sqrtf(discr);
      float invC2 = 1.0f/c2;

      float t = (-c1 - root)*invC2;
      point = rayOrigin + t*rayDir;
      E = point - coneOrigin;
      dp = dot(E, coneDir);
      if (dp > coneMinDist-MINDIST_EPS) {
        hit = true;
      }

      t = (-c1 + root)*invC2;
      point = rayOrigin + t*rayDir;
      E = point - coneOrigin;
      dp = dot(E, coneDir);
      if (dp > coneMinDist-MINDIST_EPS) {
        hit = true;
      }
    }
  }

  //  hit = true;

  if (hit)
    cudaRayPriority[idx] = 1;//(cudaDetectorRayTMax[idx] - cudaDetectorRayTMin[idx]);
  else
    cudaRayPriority[idx] = 0;

  cudaRayIds[idx] = idx;
}


__device__ void castDetectorRay(int rayIndex, int combo,
                                const float3 *cudaDetectorRayWorldOrigin,
                                const float3 *cudaDetectorRayWorldDir,
                                const float3 *cudaDetectorRayVolOrigin,
                                const float3 *cudaDetectorRayVolDir,
                                const float *cudaDetectorRayTMin,
                                const float *cudaDetectorRayTMax,
                                const float *cudaProjectionAngles,
                                int3 detectorDim,
                                int3 volDim,
                                float3 matAtten,
                                float3 matDensity,
                                float3 sourcePosition,
                                float voxelStepSize,
                                float *cudaForwardProjectionCollection) {

  float3 worldOrigin = cudaDetectorRayWorldOrigin[rayIndex];
  float3 worldDir = cudaDetectorRayWorldDir[rayIndex];
  float3 volOrigin = cudaDetectorRayVolOrigin[rayIndex];
  float3 volDir = cudaDetectorRayVolDir[rayIndex];
  float tmin = cudaDetectorRayTMin[rayIndex];
  float tmax = cudaDetectorRayTMax[rayIndex];

  if (tmin>tmax) {
    // no attenuation from the volume
    cudaForwardProjectionCollection[combo*detectorDim.x*detectorDim.y*detectorDim.z + rayIndex] = 0;
    return;
  }

  float3 volP0 = volOrigin + tmin*volDir;
  float3 volP1 = volOrigin + tmax*volDir;

  // compute ray trajectory in rotated worldspace coordinates for source falloff lookup
  // assumes detector spacing of 10 degrees!
  float theta = cudaProjectionAngles[rayIndex/(detectorDim.x*detectorDim.y)];
  float ctheta = cosf(-theta); // backwards rotation
  float stheta = sinf(-theta);
  float3 worldSourceP0 = worldOrigin + tmin*worldDir;
  float3 worldSourceP1 = worldOrigin + tmax*worldDir;
  worldSourceP0 = make_float3(ctheta*worldSourceP0.x - stheta*worldSourceP0.y,
                               stheta*worldSourceP0.x + ctheta*worldSourceP0.y,
                               worldSourceP0.z);
  worldSourceP1 = make_float3(ctheta*worldSourceP1.x - stheta*worldSourceP1.y,
                               stheta*worldSourceP1.x + ctheta*worldSourceP1.y,
                               worldSourceP1.z);
  worldSourceP0.x = -worldSourceP0.x;
  worldSourceP1.x = -worldSourceP1.x;
  worldSourceP0.z = sourcePosition.z - worldSourceP0.z;
  worldSourceP1.z = sourcePosition.z - worldSourceP1.z;

  float volLength = length(volP1-volP0);
  int numSteps = volLength / voxelStepSize + 1;


  int3 lastMat;
  float3 lastConc;
  float lastInterfaceT = tmin;

  float forwardProjection = 0;
  float sumDetectorAttenuation = 0;

  float lastInterfaceSourceAtten = sampleSourceAttenuation(combo, volP0, volDim);
  lastInterfaceSourceAtten *= tex2D(cudaSourceFalloffTexture, 
                                    worldSourceP0.x/worldSourceP0.z * (0.19/0.205) + 0.5,
                                    worldSourceP0.y/worldSourceP0.z * (0.19/0.205) + 0.5);


  for (int step=0; step<numSteps; step++) {
    float samplef = (step+0.5f) / numSteps;
    float3 volSample = lerp(volP0, volP1, samplef);

    int3 thisMat;
    float3 thisConc;
    sampleMaterialsAtPoint(combo, volSample, volDim, &thisMat, &thisConc);

    // integrate between last step and this one
    if (step > 0) {
      float interfacef;

      // interpolate an interface
      if (thisMat.x != lastMat.x) {
        interfacef = (lastConc.x-lastConc.y) / ((thisConc.x-lastConc.y) - (thisConc.y-lastConc.x));
        interfacef = clamp(interfacef, 0.0f, 1.0f);
        interfacef = lerp((step-0.5f) / numSteps, samplef, interfacef);
      }

      // fixed interface
      else {
        interfacef = (float)step / numSteps;
      }

      float interfaceT = lerp(tmin, tmax, interfacef);

      //
      // integrate this step
      //
      float tdist = interfaceT - lastInterfaceT;

      // direction to center of material segment from source
      float3 sourceRayDir = normalize((worldOrigin + ((interfaceT+lastInterfaceT)*0.5f)*worldDir) - sourcePosition);
      float ncosScatterAngle = dot(worldDir, sourceRayDir);
      float thisInterfaceSourceAtten = sampleSourceAttenuation(combo, volOrigin + interfaceT*volDir, volDim);
      
      float3 worldSourceP = lerp(worldSourceP0, worldSourceP1, interfacef);
      thisInterfaceSourceAtten *= tex2D(cudaSourceFalloffTexture,
                                        worldSourceP.x/worldSourceP.z * (0.19/0.205) + 0.5,
                                        worldSourceP.y/worldSourceP.z * (0.19/0.205) + 0.5);
      

      float density = 0;
      float massAtten = 0;
      switch (lastMat.x) {
      case 0:
        density = matDensity.x;
        massAtten = matAtten.x;
        break;
      case 1:
        density = matDensity.y;
        massAtten = matAtten.y;
        break;
      case 2:
        density = matDensity.z;
        massAtten = matAtten.z;
        break;
      }

      float voxelAttenuation = -tdist * density * massAtten;

      float attenuationFactor;
      if (voxelAttenuation == 0) {
        attenuationFactor = (lastInterfaceSourceAtten+thisInterfaceSourceAtten)*0.5f;
      }
      else {
        double a = lastInterfaceSourceAtten;
        double b = thisInterfaceSourceAtten;
        double c = voxelAttenuation;
        attenuationFactor = (exp(c) * (a+b*(c-1)) - a*(c+1)+b) / (c*c);
      }

      // attenuation between detector and start
      attenuationFactor *= exp(sumDetectorAttenuation);

      forwardProjection += (tdist *
                            sampleScatterFactor(lastMat.x, 0.5f*(1+ncosScatterAngle)) *
                            attenuationFactor);


      sumDetectorAttenuation += voxelAttenuation;

      // store info for next step
      lastInterfaceT = interfaceT;
      lastInterfaceSourceAtten = thisInterfaceSourceAtten;
    }

    lastMat = thisMat;
    lastConc = thisConc;
  }


  //
  // integrate final step
  //
  float interfaceT = tmax;
  float tdist = interfaceT - lastInterfaceT;

  // direction to center of material segment from source
  float3 sourceRayDir = normalize((worldOrigin + ((interfaceT+lastInterfaceT)*0.5f)*worldDir) - sourcePosition);
  float ncosScatterAngle = dot(worldDir, sourceRayDir);
  float thisInterfaceSourceAtten = sampleSourceAttenuation(combo, volOrigin + interfaceT*volDir, volDim);
  thisInterfaceSourceAtten *= tex2D(cudaSourceFalloffTexture,
                                    worldSourceP1.x/worldSourceP1.z * (0.19/0.205) + 0.5,
                                    worldSourceP1.y/worldSourceP1.z * (0.19/0.205) + 0.5);

  float density = 0;
  float massAtten = 0;
  switch (lastMat.x) {
  case 0:
    density = matDensity.x;
    massAtten = matAtten.x;
    break;
  case 1:
    density = matDensity.y;
    massAtten = matAtten.y;
    break;
  case 2:
    density = matDensity.z;
    massAtten = matAtten.z;
    break;
  }

  float voxelAttenuation = -tdist * density * massAtten;

  float attenuationFactor;
  if (voxelAttenuation == 0) {
    attenuationFactor = (lastInterfaceSourceAtten+thisInterfaceSourceAtten)*0.5f;
  }
  else {
    double a = lastInterfaceSourceAtten;
    double b = thisInterfaceSourceAtten;
    double c = voxelAttenuation;
    attenuationFactor = (exp(c) * (a+b*(c-1)) - a*(c+1)+b) / (c*c);
  }

  // attenuation between detector and start
  attenuationFactor *= exp(sumDetectorAttenuation);

  forwardProjection += (tdist *
                        (0.5 * (1+sampleScatterFactor(lastMat.x, ncosScatterAngle))) *
                        attenuationFactor);



  // store result
  cudaForwardProjectionCollection[combo*detectorDim.x*detectorDim.y*detectorDim.z + rayIndex] = forwardProjection;
}
                               

__global__ void castAllDetectorRays(const float3 *cudaDetectorRayWorldOrigin,
                                    const float3 *cudaDetectorRayWorldDir,
                                    const float3 *cudaDetectorRayVolOrigin,
                                    const float3 *cudaDetectorRayVolDir,
                                    const float *cudaDetectorRayTMin,
                                    const float *cudaDetectorRayTMax,
                                    const float *cudaProjectionAngles,
                                    int3 detectorDim,
                                    int3 volDim,
                                    int collectionStart,
                                    int collectionEnd,
                                    float3 matAtten,
                                    float3 matDensity,
                                    float3 sourcePosition,
                                    float voxelStepSize,
                                    float *cudaForwardProjectionCollection) {

  unsigned int gx = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int gy = blockIdx.y*blockDim.y + threadIdx.y;
  unsigned int gz = blockIdx.z*blockDim.z + threadIdx.z;
  if (gx>=detectorDim.x || gy>=detectorDim.y || 
      gz<detectorDim.z*collectionStart || gz>=detectorDim.z*collectionEnd)
    return;
  unsigned int rayIndex = (gz%detectorDim.z)*detectorDim.x*detectorDim.y + gy*detectorDim.x + gx;
  int combo = gz/detectorDim.z;

  castDetectorRay(rayIndex, combo,
                  cudaDetectorRayWorldOrigin,
                  cudaDetectorRayWorldDir,
                  cudaDetectorRayVolOrigin,
                  cudaDetectorRayVolDir,
                  cudaDetectorRayTMin,
                  cudaDetectorRayTMax,
                  cudaProjectionAngles,
                  detectorDim,
                  volDim,
                  matAtten,
                  matDensity,
                  sourcePosition,
                  voxelStepSize,
                  cudaForwardProjectionCollection);
}



__global__ void castPrioritizedDetectorRays(const unsigned int *cudaRayIds,
                                            const float *cudaRayPriority,
                                            const float *cudaCurrentForwardProjection,
                                            const float3 *cudaDetectorRayWorldOrigin,
                                            const float3 *cudaDetectorRayWorldDir,
                                            const float3 *cudaDetectorRayVolOrigin,
                                            const float3 *cudaDetectorRayVolDir,
                                            const float *cudaDetectorRayTMin,
                                            const float *cudaDetectorRayTMax,
                                            const float *cudaProjectionAngles,
                                            int3 detectorDim,
                                            int3 volDim,
                                            int collectionStart,
                                            int collectionEnd,
                                            float3 matAtten,
                                            float3 matDensity,
                                            float3 sourcePosition,
                                            float voxelStepSize,
                                            float *cudaForwardProjectionCollection) {

  unsigned int gx = blockIdx.x*blockDim.x + threadIdx.x;
  if (gx<detectorDim.x*detectorDim.y*detectorDim.z*collectionStart ||
      gx>=detectorDim.x*detectorDim.y*detectorDim.z*collectionEnd)
    return;

  //unsigned int rayIndexIndex = gx % (detectorDim.x*detectorDim.y*detectorDim.z);
  //int combo = gx / (detectorDim.x*detectorDim.y*detectorDim.z);
  unsigned int rayIndexIndex = gx % (detectorDim.x*detectorDim.y*detectorDim.z);
  int combo = gx / (detectorDim.x*detectorDim.y*detectorDim.z);

  unsigned int rayIndex = cudaRayIds[rayIndexIndex];
  float priority = cudaRayPriority[rayIndexIndex];

  if (priority == 0) {
    cudaForwardProjectionCollection[combo*
                                    detectorDim.x*
                                    detectorDim.y*
                                    detectorDim.z +
                                    rayIndex] = cudaCurrentForwardProjection[rayIndex];
  }
  else {
    castDetectorRay(rayIndex, combo,
                    cudaDetectorRayWorldOrigin,
                    cudaDetectorRayWorldDir,
                    cudaDetectorRayVolOrigin,
                    cudaDetectorRayVolDir,
                    cudaDetectorRayTMin,
                    cudaDetectorRayTMax,
                    cudaProjectionAngles,
                    detectorDim,
                    volDim,
                    matAtten,
                    matDensity,
                    sourcePosition,
                    voxelStepSize,
                    cudaForwardProjectionCollection);
  }
}




void MarkovContext::CudaForwardProject(int collectionSize) const {
  TIMER_START("CudaForwardProject()");

  int3 volDim = make_int3(mGeometry.GetVolumeNodeSamplesX(),
                          mGeometry.GetVolumeNodeSamplesY(),
                          mGeometry.GetVolumeNodeSamplesZ());
  int3 detectorDim = make_int3(mGeometry.GetDetectorSamplesWidth(),
                               mGeometry.GetDetectorSamplesHeight(),
                               mGeometry.GetNumProjectionAngles());

  float3 matAtten = make_float3(mMaterials[0].GetMassAttenuationCoefficient(),
                                mMaterials[1].GetMassAttenuationCoefficient(),
                                mMaterials[2].GetMassAttenuationCoefficient());
  float3 matDensity = make_float3(mMaterials[0].GetDensity(),
                                  mMaterials[1].GetDensity(),
                                  mMaterials[2].GetDensity());

  float3 sourcePosition = make_float3(mGeometry.GetSourcePosition()[0],
                                      mGeometry.GetSourcePosition()[1],
                                      mGeometry.GetSourcePosition()[2]);
                                

  dim3 dimBlock(8, 8, 1);
  dim3 dimGrid(1+(detectorDim.x-1) / dimBlock.x, 
               1+(detectorDim.y-1) / dimBlock.y, 
               1+(collectionSize*detectorDim.z-1) / dimBlock.z);

  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, collectionSize, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    castAllDetectorRays<<<dimGrid, dimBlock>>>(cudaDetectorRayWorldOrigin[g],
                                               cudaDetectorRayWorldDir[g],
                                               cudaDetectorRayVolOrigin[g],
                                               cudaDetectorRayVolDir[g],
                                               cudaDetectorRayTMin[g],
                                               cudaDetectorRayTMax[g],
                                               cudaProjectionAngles[g],
                                               detectorDim,
                                               volDim,
                                               collectionStart,
                                               collectionEnd,
                                               matAtten,
                                               matDensity,
                                               sourcePosition,
                                               mVoxelStepSize,
                                               cudaForwardProjectionCollection[g]);
  }

  TIMER_STOP;
}


#ifdef CUDA_ENABLE_UPDATE_FORWARD_PROJECTION
void MarkovContext::CudaUpdateForwardProjection(int collectionSize,
                                                const Cone &attenChangeCone) const {

  if (mGpuIds.size() > 1) {
    // currently, forward projections are not copied between gpus when a config is accepted
    std::cerr<<"CudaUpdateForwardProjection() currently does not support more than one GPU!"<<std::endl;
    return;
  }

  TIMER_START("CudaUpdateForwardProjection()");

  int3 volDim = make_int3(mGeometry.GetVolumeNodeSamplesX(),
                          mGeometry.GetVolumeNodeSamplesY(),
                          mGeometry.GetVolumeNodeSamplesZ());
  int3 detectorDim = make_int3(mGeometry.GetDetectorSamplesWidth(),
                               mGeometry.GetDetectorSamplesHeight(),
                               mGeometry.GetNumProjectionAngles());

  float3 matAtten = make_float3(mMaterials[0].GetMassAttenuationCoefficient(),
                                mMaterials[1].GetMassAttenuationCoefficient(),
                                mMaterials[2].GetMassAttenuationCoefficient());
  float3 matDensity = make_float3(mMaterials[0].GetDensity(),
                                  mMaterials[1].GetDensity(),
                                  mMaterials[2].GetDensity());

  float3 sourcePosition = make_float3(mGeometry.GetSourcePosition()[0],
                                      mGeometry.GetSourcePosition()[1],
                                      mGeometry.GetSourcePosition()[2]);

  // prioritize each ray
  dim3 dimBlockPrioritize(32,1,1);
  dim3 dimGridPrioritize(1+(mGeometry.GetTotalProjectionSamples()-1) / dimBlockPrioritize.x, 1, 1);
  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    prioritizeDetectorRays<<<dimGridPrioritize, dimBlockPrioritize>>>(cudaDetectorRayWorldOrigin[g],
                                                                      cudaDetectorRayWorldDir[g],
                                                                      cudaDetectorRayTMin[g],
                                                                      cudaDetectorRayTMax[g],
                                                                      mGeometry.GetTotalProjectionSamples(),
                                                                      make_float3(attenChangeCone.mOrigin[0],
                                                                                  attenChangeCone.mOrigin[1],
                                                                                  attenChangeCone.mOrigin[2]),
                                                                      make_float3(attenChangeCone.mDir[0],
                                                                                  attenChangeCone.mDir[1],
                                                                                  attenChangeCone.mDir[2]),
                                                                      attenChangeCone.mCosTheta,
                                                                      attenChangeCone.mMinDist,
                                                                      cudaRayIds[g],
                                                                      cudaRayPriority[g]);
  }

  // sort rays by priority
  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    thrust::device_ptr<unsigned int> thrustIds = thrust::device_pointer_cast(cudaRayIds[g]);
    thrust::device_ptr<float> thrustPriorities = thrust::device_pointer_cast(cudaRayPriority[g]);
    thrust::sort_by_key(thrustPriorities, thrustPriorities+mGeometry.GetTotalProjectionSamples(), thrustIds);
    //thrust::stable_sort_by_key(thrustPriorities, thrustPriorities+mGeometry.GetTotalProjectionSamples(), thrustIds);

    /*
      vector<float> priorities(mGeometry.GetTotalProjectionSamples());
      checkCudaErrors(cudaMemcpy(&priorities[0], cudaRayPriority, sizeof(float)*mGeometry.GetTotalProjectionSamples(), cudaMemcpyDeviceToHost));

      int skippedRays = 0;
      for (int i=0; i<priorities.size(); i++) {
      if (priorities[i] == 0) {
      skippedRays++;
      }
      }

      std::cerr<<"skipped "<<skippedRays<<" of "<<priorities.size()<<std::endl;
    */


    // cast only non-zero priority rays
    /*
      dim3 dimBlock(8, 8, 1);
      dim3 dimGrid(1+(detectorDim.x-1) / dimBlock.x, 
      1+(detectorDim.y-1) / dimBlock.y, 
      1+(collectionSize*detectorDim.z-1) / dimBlock.z);
    */
  }


  dim3 dimBlock(32,1,1);
  dim3 dimGrid(1+(collectionSize * mGeometry.GetTotalProjectionSamples()-1) / dimBlockPrioritize.x, 1, 1);

  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, collectionSize, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    castPrioritizedDetectorRays<<<dimGrid, dimBlock>>>(cudaRayIds[g],
                                                       cudaRayPriority[g],
                                                       cudaCurrentForwardProjection[g],
                                                       cudaDetectorRayWorldOrigin[g],
                                                       cudaDetectorRayWorldDir[g],
                                                       cudaDetectorRayVolOrigin[g],
                                                       cudaDetectorRayVolDir[g],
                                                       cudaDetectorRayTMin[g],
                                                       cudaDetectorRayTMax[g],
                                                       cudaProjectionAngles[g],
                                                       detectorDim,
                                                       volDim,
                                                       collectionStart,
                                                       collectionEnd,
                                                       matAtten,
                                                       matDensity,
                                                       sourcePosition,
                                                       mVoxelStepSize,
                                                       cudaForwardProjectionCollection[g]);
  }

  TIMER_STOP;
}
#endif

void MarkovContext::CudaGetForwardProjection(int collectionSize,
                                             vector< vector<float> > &forwardProjectionCollection) const {
  TIMER_START("CudaGetForwardProjection()");

  // copy results back
  //vector<float> outputData(collectionSize * mGeometry.GetTotalProjectionSamples());
  float *outputData;
  checkCudaErrors(cudaHostAlloc(&outputData,
                                collectionSize * mGeometry.GetTotalProjectionSamples() * sizeof(float),
                                cudaHostAllocDefault));

  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, collectionSize, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    if ((collectionEnd-collectionStart) > 0) {
      checkCudaErrors(cudaMemcpyAsync(&outputData[collectionStart * mGeometry.GetTotalProjectionSamples()],
                                      &cudaForwardProjectionCollection[g][collectionStart * mGeometry.GetTotalProjectionSamples()],
                                      (collectionEnd-collectionStart) * mGeometry.GetTotalProjectionSamples()*sizeof(float),
                                      cudaMemcpyDeviceToHost));
    }
  }

  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    checkCudaErrors(cudaDeviceSynchronize());
  }


  forwardProjectionCollection.resize(collectionSize);
  for (int c=0; c<collectionSize; c++) {
    forwardProjectionCollection[c].resize(mGeometry.GetTotalProjectionSamples());
    for (int i=0; i<mGeometry.GetTotalProjectionSamples(); i++) {
      forwardProjectionCollection[c][i] = outputData[c*mGeometry.GetTotalProjectionSamples() + i];
    }
  }

  checkCudaErrors(cudaFreeHost(outputData));


  TIMER_STOP;
}



//==================================================================================================
//==================================================================================================
//==================================================================================================
__global__ void projectionToError(int totalProjectionSamples,
                                  int collectionSize,
                                  const float *baselineProjection,
                                  const float *forwardProjection,
                                  float *forwardProjectionError) {

  unsigned int cpi = blockIdx.x*blockDim.x + threadIdx.x;
  if (cpi >= totalProjectionSamples*collectionSize)
    return;

  int pi = cpi % totalProjectionSamples;
  float df = forwardProjection[cpi] - baselineProjection[pi];
  forwardProjectionError[cpi] = df*df;
}


void MarkovContext::CudaGetProjectionError(int collectionSize, vector<float> &errors) const {

  TIMER_START("CudaGetProjectionError()");

  // compute squared errors
  int totalProjectionSamples = mGeometry.GetTotalProjectionSamples();
  dim3 dimBlock(32, 1, 1);
  dim3 dimGrid(1+((totalProjectionSamples*collectionSize)-1) / dimBlock.x, 1, 1);

  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    projectionToError<<<dimGrid, dimBlock>>>(totalProjectionSamples,
                                             collectionSize,
                                             cudaBaselineProjection[g],
                                             cudaForwardProjectionCollection[g],
                                             cudaForwardProjectionError[g]);
  }

  // use thrust to sum the errors for each material combo
  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    thrust::device_ptr<float> dev_ptr = thrust::device_pointer_cast(cudaForwardProjectionError[g]);
    errors.resize(collectionSize);
    for (int c=0; c<collectionSize; c++) {
      errors[c] = (thrust::reduce(dev_ptr+c*totalProjectionSamples,
                                  dev_ptr+(c+1)*totalProjectionSamples)
                   * mGeometry.GetDetectorPixelArea());
    }
  }

  TIMER_STOP;
}



//==================================================================================================
//==================================================================================================
//==================================================================================================
__global__ void cudaUpdateVolume1(float4 *vol,
                                  int idx,
                                  int mat) {

  // only a single thread needs to do anything
  if (blockIdx.x != 0 || threadIdx.x != 0)
    return;

  float4 nv = make_float4(0,0,0,0);
  switch (mat) {
  case 0:  nv.x = 1;  break;
  case 1:  nv.y = 1;  break;
  case 2:  nv.z = 1;  break;
  }

  vol[idx] = nv;
}


__global__ void cudaUpdateVolume2(float4 *vol,
                                  int idx,
                                  int mat,
                                  int idx2,
                                  int mat2) {

  // only a single thread needs to do anything
  if (blockIdx.x != 0 || threadIdx.x != 0)
    return;

  float4 nv = make_float4(0,0,0,0);
  switch (mat) {
  case 0:  nv.x = 1;  break;
  case 1:  nv.y = 1;  break;
  case 2:  nv.z = 1;  break;
  }
  vol[idx] = nv;

  float4 nv2 = make_float4(0,0,0,0);
  switch (mat2) {
  case 0:  nv2.x = 1;  break;
  case 1:  nv2.y = 1;  break;
  case 2:  nv2.z = 1;  break;
  }
  vol[idx2] = nv2;
}



void MarkovContext::CudaAcceptNextConfig(const GibbsProposal &proposal, int c) const {
  TIMER_START("CudaAcceptNextConfig()");

  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, NUM_MATERIALS*NUM_MATERIALS, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));

    // update all collection volumes
    for (int c2=collectionStart; c2<collectionEnd; c2++) {
      float4 *cudaVolume = cudaVolumeLinearCollection[g] + c2*mGeometry.GetTotalVolumeNodeSamples();

      // apply changes
      // first proposal only
      if (proposal.first>=0 && proposal.second<0) {
        dim3 dimBlock(1,1,1);
        dim3 dimGrid(1,1,1);
        cudaUpdateVolume1<<<dimGrid, dimBlock>>>(cudaVolume,
                                                 proposal.first,
                                                 c%NUM_MATERIALS);
      }

      // apply changes
      // both proposals
      else if (proposal.first>=0 && proposal.second>=0) {
        dim3 dimBlock(1,1,1);
        dim3 dimGrid(1,1,1);
        cudaUpdateVolume2<<<dimGrid, dimBlock>>>(cudaVolume,
                                                 proposal.first,
                                                 c%NUM_MATERIALS,
                                                 proposal.second,
                                                 c/NUM_MATERIALS);
      }
    }


    // set the accepted forward projection as current
#ifdef CUDA_ENABLE_UPDATE_FORWARD_PROJECTION
    checkCudaErrors(cudaMemcpyAsync(cudaCurrentForwardProjection[g],
                                    cudaForwardProjectionCollection[g] + c*mGeometry.GetTotalProjectionSamples(),
                                    sizeof(float)*mGeometry.GetTotalProjectionSamples(),
                                    cudaMemcpyDeviceToDevice));
#endif
  }

  TIMER_STOP;
}


void MarkovContext::CudaSetCurrentVolume(const vector<unsigned char> &matids) const {

  // setup individual channels per material
  vector<float4> volumeData(mGeometry.GetTotalVolumeNodeSamples());
  memset(&volumeData[0], 0, sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples());
  for (int i=0; i<mGeometry.GetTotalVolumeNodeSamples(); i++) {
    switch (matids[i]) {
    case 0:  volumeData[i].x = 1;  break;
    case 1:  volumeData[i].y = 1;  break;
    case 2:  volumeData[i].z = 1;  break;
    }
  }

  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));

    // set all collection volumes
    for (int c=0; c<NUM_MATERIALS*NUM_MATERIALS; c++) {
      checkCudaErrors(cudaMemcpy(cudaVolumeLinearCollection[g] + c*mGeometry.GetTotalVolumeNodeSamples(),
                                 &volumeData[0],
                                 sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples(),
                                 cudaMemcpyHostToDevice));
    }
  }
}


void MarkovContext::CudaGetVolumeCollection(vector< vector<unsigned char> > &volumeCollection) const {
  volumeCollection.resize(NUM_MATERIALS*NUM_MATERIALS);
  for (int c=0; c<NUM_MATERIALS*NUM_MATERIALS; c++) {

    vector<float4> fvol(mGeometry.GetTotalVolumeNodeSamples());

    checkCudaErrors(cudaMemcpy(&fvol[0],
                               cudaVolumeLinearCollection + c*mGeometry.GetTotalVolumeNodeSamples(),
                               sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples(),
                               cudaMemcpyDeviceToHost));


    volumeCollection[c].resize(mGeometry.GetTotalVolumeNodeSamples());
    for (int i=0; i<mGeometry.GetTotalVolumeNodeSamples(); i++) {
      if (fvol[i].x == 1)
        volumeCollection[c][i] = 0;
      else if (fvol[i].y == 1)
        volumeCollection[c][i] = 1;
      else if (fvol[i].z == 1)
        volumeCollection[c][i] = 2;
      else {
        std::cerr<<"bogus volume!"<<std::endl;
        exit(0);
      }
    }
    
  }
}


void MarkovContext::CudaSetVolumeCollection(const GibbsProposal &proposal) const {
  TIMER_START("CudaSetVolumeCollection()");

  for (int g=0; g<mGpuIds.size(); g++) {
    int collectionStart, collectionEnd;
    CudaGetCollectionStartEnd(g, NUM_MATERIALS*NUM_MATERIALS, collectionStart, collectionEnd);

    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
 
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
    cudaExtent volumeSize = make_cudaExtent(mGeometry.GetVolumeNodeSamplesX(),
                                            mGeometry.GetVolumeNodeSamplesY(),
                                            mGeometry.GetVolumeNodeSamplesZ());

    // apply changes to each collection volume
    for (int c=collectionStart; c<collectionEnd; c++) {
      float4 *cudaVolume = cudaVolumeLinearCollection[g] + c*mGeometry.GetTotalVolumeNodeSamples();

      // apply changes
      // first proposal only
      if (proposal.first>=0 && proposal.second<0) {
        dim3 dimBlock(1,1,1);
        dim3 dimGrid(1,1,1);
        cudaUpdateVolume1<<<dimGrid, dimBlock>>>(cudaVolume,
                                                 proposal.first,
                                                 c%NUM_MATERIALS);
      }

      // apply changes
      // both proposals
      else if (proposal.first>=0 && proposal.second>=0) {
        dim3 dimBlock(1,1,1);
        dim3 dimGrid(1,1,1);
        cudaUpdateVolume2<<<dimGrid, dimBlock>>>(cudaVolume,
                                                 proposal.first,
                                                 c%NUM_MATERIALS,
                                                 proposal.second,
                                                 c/NUM_MATERIALS);
      }
    }


    // copy linear collection volumes to arrays
    for (int c=collectionStart; c<collectionEnd; c++) {
      float4 *cudaVolume = cudaVolumeLinearCollection[g] + c*mGeometry.GetTotalVolumeNodeSamples();

      cudaArray **cudaArray = NULL;
      switch (c) {
      case 0:
        cudaArray = &cudaVolumeArray00[g];
        break;
      case 1:
        cudaArray = &cudaVolumeArray01[g];
        break;
      case 2:
        cudaArray = &cudaVolumeArray02[g];
        break;
      case 3:
        cudaArray = &cudaVolumeArray10[g];
        break;
      case 4:
        cudaArray = &cudaVolumeArray11[g];
        break;
      case 5:
        cudaArray = &cudaVolumeArray12[g];
        break;
      case 6:
        cudaArray = &cudaVolumeArray20[g];
        break;
      case 7:
        cudaArray = &cudaVolumeArray21[g];
        break;
      case 8:
        cudaArray = &cudaVolumeArray22[g];
        break;
      }

      // copy data to 3D array
      cudaMemcpy3DParms copyParams = {0};
      copyParams.srcPtr   = make_cudaPitchedPtr(cudaVolume, volumeSize.width*sizeof(float4), volumeSize.width, volumeSize.height);
      copyParams.dstArray = *cudaArray;
      copyParams.extent   = volumeSize;
      copyParams.kind     = cudaMemcpyDeviceToDevice;
      checkCudaErrors(cudaMemcpy3DAsync(&copyParams));
    }
  }

  TIMER_STOP;
}


//==================================================================================================
//==================================================================================================
//==================================================================================================
void MarkovContext::CudaSetBaselineProjection() const {
  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));

    checkCudaErrors(cudaMalloc(&cudaBaselineProjection[g], sizeof(float)*mBaselineProjection.size()));
    checkCudaErrors(cudaMemcpy(cudaBaselineProjection[g], &mBaselineProjection[0], sizeof(float)*mBaselineProjection.size(), cudaMemcpyHostToDevice));
  }
}


template <typename T>
void SetTextureParams(T *cudaTexture) {
  cudaTexture->normalized = true;
  cudaTexture->filterMode = cudaFilterModeLinear;
  cudaTexture->addressMode[0] = cudaAddressModeClamp;
  cudaTexture->addressMode[1] = cudaAddressModeClamp;
  cudaTexture->addressMode[2] = cudaAddressModeClamp;
}

void MarkovContext::CudaInitialize() {

  mGpuIds.clear();

  int gpu_bitfield = mGeometry.GetGPUBitfield();
  for (int g=0; g<MAX_GPUS; g++) {
    if ((1<<g)&gpu_bitfield) {
      int id = gpuDeviceInit(g);
      if (id>=0) {
        mGpuIds.push_back(id);

        // get device name
        cudaDeviceProp deviceProps;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProps, id));
        printf("CUDA device [%s]\n", deviceProps.name);
      }
    }
  }

  if (mGpuIds.empty()) {
    std::cerr<<"No GPUs selected!!"<<std::endl;
  }


  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));

    // initialize material data
    for (int m=0; m<mMaterials.size(); m++) {
      vector<float> scatterFactors;
      mMaterials[m].GetScatterFactorArray(scatterFactors);

      // Allocate array and copy image data
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
      cudaArray *cuArray;
      checkCudaErrors(cudaMallocArray(&cuArray,
                                      &channelDesc,
                                      (int)scatterFactors.size(),
                                      1,
                                      cudaArrayDefault));
      cudaMaterialArray[g][m] = cuArray;

      checkCudaErrors(cudaMemcpyToArray(cuArray,
                                        0,
                                        0,
                                        &scatterFactors[0],
                                        (int)scatterFactors.size() * sizeof(float),
                                        cudaMemcpyHostToDevice));


      texture<float, 1, cudaReadModeElementType> *cudaMaterialTexture = NULL;
      switch (m) {
      case 0:
        cudaMaterialTexture = &cudaMaterialTextures0; break;
      case 1:
        cudaMaterialTexture = &cudaMaterialTextures1; break;
      case 2:
        cudaMaterialTexture = &cudaMaterialTextures2; break;
      }

      // Set texture parameters
      cudaMaterialTexture->addressMode[0] = cudaAddressModeClamp;
      cudaMaterialTexture->addressMode[1] = cudaAddressModeClamp;
      cudaMaterialTexture->filterMode = cudaFilterModeLinear;
      cudaMaterialTexture->normalized = true;    // access with normalized texture coordinates
    
      // Bind the array to the texture
      checkCudaErrors(cudaBindTextureToArray(*cudaMaterialTexture, cuArray, channelDesc));
    }


    // volume data
    cudaExtent volumeSize = make_cudaExtent(mGeometry.GetVolumeNodeSamplesX(),
                                            mGeometry.GetVolumeNodeSamplesY(),
                                            mGeometry.GetVolumeNodeSamplesZ());
    checkCudaErrors(cudaMalloc(&cudaVolumeLinearCollection[g], sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples()*NUM_MATERIALS*NUM_MATERIALS));
    cudaChannelFormatDesc channelDesc4 = cudaCreateChannelDesc<float4>();
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray00[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray01[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray02[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray10[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray11[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray12[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray20[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray21[g], &channelDesc4, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray22[g], &channelDesc4, volumeSize));

    SetTextureParams(&cudaVolumeTextures00);
    SetTextureParams(&cudaVolumeTextures01);
    SetTextureParams(&cudaVolumeTextures02);
    SetTextureParams(&cudaVolumeTextures10);
    SetTextureParams(&cudaVolumeTextures11);
    SetTextureParams(&cudaVolumeTextures12);
    SetTextureParams(&cudaVolumeTextures20);
    SetTextureParams(&cudaVolumeTextures21);
    SetTextureParams(&cudaVolumeTextures22);
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures00, cudaVolumeArray00[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures01, cudaVolumeArray01[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures02, cudaVolumeArray02[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures10, cudaVolumeArray10[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures11, cudaVolumeArray11[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures12, cudaVolumeArray12[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures20, cudaVolumeArray20[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures21, cudaVolumeArray21[g], channelDesc4));
    checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures22, cudaVolumeArray22[g], channelDesc4));


    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray00[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray01[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray02[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray10[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray11[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray12[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray20[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray21[g], &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray22[g], &channelDesc, volumeSize));

    SetTextureParams(&cudaSourceTextures00);
    SetTextureParams(&cudaSourceTextures01);
    SetTextureParams(&cudaSourceTextures02);
    SetTextureParams(&cudaSourceTextures10);
    SetTextureParams(&cudaSourceTextures11);
    SetTextureParams(&cudaSourceTextures12);
    SetTextureParams(&cudaSourceTextures20);
    SetTextureParams(&cudaSourceTextures21);
    SetTextureParams(&cudaSourceTextures22);
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures00, cudaSourceArray00[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures01, cudaSourceArray01[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures02, cudaSourceArray02[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures10, cudaSourceArray10[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures11, cudaSourceArray11[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures12, cudaSourceArray12[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures20, cudaSourceArray20[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures21, cudaSourceArray21[g], channelDesc));
    checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures22, cudaSourceArray22[g], channelDesc));


    // upload all of the ray info
    checkCudaErrors(cudaMalloc(&cudaSourceRayVolOrigin[g], sizeof(float3)*mSourceRayVolOrigin.size()));
    checkCudaErrors(cudaMemcpy(cudaSourceRayVolOrigin[g], &mSourceRayVolOrigin[0], sizeof(float3)*mSourceRayVolOrigin.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaSourceRayVolDir[g], sizeof(float3)*mSourceRayVolDir.size()));
    checkCudaErrors(cudaMemcpy(cudaSourceRayVolDir[g], &mSourceRayVolDir[0], sizeof(float3)*mSourceRayVolDir.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaSourceRayTMin[g], sizeof(float)*mSourceRayTMin.size()));
    checkCudaErrors(cudaMemcpy(cudaSourceRayTMin[g], &mSourceRayTMin[0], sizeof(float)*mSourceRayTMin.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaSourceRayTMax[g], sizeof(float)*mSourceRayTMax.size()));
    checkCudaErrors(cudaMemcpy(cudaSourceRayTMax[g], &mSourceRayTMax[0], sizeof(float)*mSourceRayTMax.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaDetectorRayVolOrigin[g], sizeof(float3)*mDetectorRayVolOrigin.size()));
    checkCudaErrors(cudaMemcpy(cudaDetectorRayVolOrigin[g], &mDetectorRayVolOrigin[0], sizeof(float3)*mDetectorRayVolOrigin.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaDetectorRayVolDir[g], sizeof(float3)*mDetectorRayVolDir.size()));
    checkCudaErrors(cudaMemcpy(cudaDetectorRayVolDir[g], &mDetectorRayVolDir[0], sizeof(float3)*mDetectorRayVolDir.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaDetectorRayWorldOrigin[g], sizeof(float3)*mDetectorRayOrigin.size()));
    checkCudaErrors(cudaMemcpy(cudaDetectorRayWorldOrigin[g], &mDetectorRayOrigin[0], sizeof(float3)*mDetectorRayOrigin.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaDetectorRayWorldDir[g], sizeof(float3)*mDetectorRayDir.size()));
    checkCudaErrors(cudaMemcpy(cudaDetectorRayWorldDir[g], &mDetectorRayDir[0], sizeof(float3)*mDetectorRayDir.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaDetectorRayTMin[g], sizeof(float)*mDetectorRayTMin.size()));
    checkCudaErrors(cudaMemcpy(cudaDetectorRayTMin[g], &mDetectorRayTMin[0], sizeof(float)*mDetectorRayTMin.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc(&cudaDetectorRayTMax[g], sizeof(float)*mDetectorRayTMax.size()));
    checkCudaErrors(cudaMemcpy(cudaDetectorRayTMax[g], &mDetectorRayTMax[0], sizeof(float)*mDetectorRayTMax.size(), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc(&cudaRayIds[g], sizeof(unsigned int)*mGeometry.GetTotalProjectionSamples()));
    checkCudaErrors(cudaMalloc(&cudaRayPriority[g], sizeof(float)*mGeometry.GetTotalProjectionSamples()));
    checkCudaErrors(cudaMalloc(&cudaCurrentForwardProjection[g], sizeof(float)*mGeometry.GetTotalProjectionSamples()));

    vector<float> projectionAngles;
    for (int i=0; i<mGeometry.GetNumProjectionAngles(); i++) {
      projectionAngles.push_back(mGeometry.GetProjectionAngle(i));
    }
    checkCudaErrors(cudaMalloc(&cudaProjectionAngles[g], sizeof(float)*mGeometry.GetNumProjectionAngles()));
    checkCudaErrors(cudaMemcpy(cudaProjectionAngles[g], &projectionAngles[0], sizeof(float)*mGeometry.GetNumProjectionAngles(), cudaMemcpyHostToDevice));


    // precompute some source attenuation info
    vector<float> sourceScale(mGeometry.GetTotalVolumeNodeSamples());
    for (int nvi=0; nvi<mGeometry.GetTotalVolumeNodeSamples(); nvi++) {
      int x,y,z;
      mGeometry.VolumeIndexToNodeCoord(nvi, x,y,z);

      Vec3f voxelPosition;
      mGeometry.VolumeToWorld(Vec3f((float)x,(float)y,(float)z), voxelPosition);

      Vec3f diff = voxelPosition - mGeometry.GetSourcePosition();
      float maxt = diff.Length();
      Vec3f dir = diff / maxt;
    
      sourceScale[nvi] = 1 / (maxt*maxt);
    }

    checkCudaErrors(cudaMalloc(&cudaSourceScale[g], sizeof(float)*sourceScale.size()));
    checkCudaErrors(cudaMemcpy(cudaSourceScale[g], &sourceScale[0], sizeof(float)*sourceScale.size(), cudaMemcpyHostToDevice));


    // source falloff texture
    cudaChannelFormatDesc falloffChannelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    checkCudaErrors(cudaMallocArray(&cudaFalloffArray[g],
                                    &falloffChannelDesc,
                                    mGeometry.GetSourceAttenMapWidth(),
                                    mGeometry.GetSourceAttenMapHeight(),
                                    cudaArrayDefault));
    checkCudaErrors(cudaMemcpyToArray(cudaFalloffArray[g],
                                      0,
                                      0,
                                      &mGeometry.GetSourceAttenMap()[0],
                                      mGeometry.GetSourceAttenMapWidth()*mGeometry.GetSourceAttenMapHeight() * sizeof(float),
                                      cudaMemcpyHostToDevice));
    cudaSourceFalloffTexture.addressMode[0] = cudaAddressModeClamp;
    cudaSourceFalloffTexture.addressMode[1] = cudaAddressModeClamp;
    cudaSourceFalloffTexture.filterMode = cudaFilterModeLinear;
    cudaSourceFalloffTexture.normalized = true;    // access with normalized texture coordinates
    checkCudaErrors(cudaBindTextureToArray(cudaSourceFalloffTexture, cudaFalloffArray[g], falloffChannelDesc));



    checkCudaErrors(cudaMalloc(&cudaSourceAttenuationCollection[g], sizeof(float)*mGeometry.GetTotalVolumeNodeSamples() * NUM_MATERIALS*NUM_MATERIALS));
    checkCudaErrors(cudaMalloc(&cudaForwardProjectionCollection[g], sizeof(float)*mGeometry.GetTotalProjectionSamples() * NUM_MATERIALS*NUM_MATERIALS));
    checkCudaErrors(cudaMalloc(&cudaForwardProjectionError[g], sizeof(float)*mGeometry.GetTotalProjectionSamples() * NUM_MATERIALS*NUM_MATERIALS));
  }
}


void MarkovContext::CudaShutdown() const {
  for (int g=0; g<mGpuIds.size(); g++) {
    checkCudaErrors(cudaSetDevice(mGpuIds[g]));
    
    checkCudaErrors(cudaFree(cudaBaselineProjection[g]));  cudaBaselineProjection[g]=NULL;

    for (int m=0; m<mMaterials.size(); m++) {
      checkCudaErrors(cudaFree(cudaMaterialArray[g][m]));  cudaMaterialArray[g][m]=NULL;
    }

    checkCudaErrors(cudaFree(cudaVolumeLinearCollection[g]));  cudaVolumeLinearCollection[g]=NULL;

    checkCudaErrors(cudaFreeArray(cudaVolumeArray00[g]));  cudaVolumeArray00[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray01[g]));  cudaVolumeArray01[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray02[g]));  cudaVolumeArray02[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray10[g]));  cudaVolumeArray10[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray11[g]));  cudaVolumeArray11[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray12[g]));  cudaVolumeArray12[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray20[g]));  cudaVolumeArray20[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray21[g]));  cudaVolumeArray21[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaVolumeArray22[g]));  cudaVolumeArray22[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray00[g]));  cudaSourceArray00[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray01[g]));  cudaSourceArray01[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray02[g]));  cudaSourceArray02[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray10[g]));  cudaSourceArray10[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray11[g]));  cudaSourceArray11[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray12[g]));  cudaSourceArray12[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray20[g]));  cudaSourceArray20[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray21[g]));  cudaSourceArray21[g]=NULL;
    checkCudaErrors(cudaFreeArray(cudaSourceArray22[g]));  cudaSourceArray22[g]=NULL;

    checkCudaErrors(cudaFree(cudaSourceRayVolOrigin[g]));  cudaSourceRayVolOrigin[g]=NULL;
    checkCudaErrors(cudaFree(cudaSourceRayVolDir[g]));  cudaSourceRayVolDir[g]=NULL;
    checkCudaErrors(cudaFree(cudaSourceRayTMin[g]));  cudaSourceRayTMin[g]=NULL;
    checkCudaErrors(cudaFree(cudaSourceRayTMax[g]));  cudaSourceRayTMax[g]=NULL;
    checkCudaErrors(cudaFree(cudaDetectorRayVolOrigin[g]));  cudaDetectorRayVolOrigin[g]=NULL;
    checkCudaErrors(cudaFree(cudaDetectorRayVolDir[g]));  cudaDetectorRayVolDir[g]=NULL;
    checkCudaErrors(cudaFree(cudaDetectorRayWorldOrigin[g]));  cudaDetectorRayWorldOrigin[g]=NULL;
    checkCudaErrors(cudaFree(cudaDetectorRayWorldDir[g]));  cudaDetectorRayWorldDir[g]=NULL;
    checkCudaErrors(cudaFree(cudaDetectorRayTMin[g]));  cudaDetectorRayTMin[g]=NULL;
    checkCudaErrors(cudaFree(cudaDetectorRayTMax[g]));  cudaDetectorRayTMax[g]=NULL;

    checkCudaErrors(cudaFree(cudaRayIds[g]));  cudaRayIds[g]=NULL;
    checkCudaErrors(cudaFree(cudaRayPriority[g]));  cudaRayPriority[g]=NULL;
    checkCudaErrors(cudaFree(cudaCurrentForwardProjection[g]));  cudaCurrentForwardProjection[g]=NULL;

    checkCudaErrors(cudaFree(cudaProjectionAngles[g]));  cudaProjectionAngles[g]=NULL;

    checkCudaErrors(cudaFree(cudaSourceScale[g]));  cudaSourceScale[g]=NULL;

    checkCudaErrors(cudaFreeArray(cudaFalloffArray[g]));  cudaFalloffArray[g]=NULL;

    checkCudaErrors(cudaFree(cudaSourceAttenuationCollection[g]));  cudaSourceAttenuationCollection[g]=NULL;
    checkCudaErrors(cudaFree(cudaForwardProjectionCollection[g]));  cudaForwardProjectionCollection[g]=NULL;
    checkCudaErrors(cudaFree(cudaForwardProjectionError[g]));  cudaForwardProjectionError[g]=NULL;

    cudaDeviceReset();
  }
}


void MarkovContext::CudaGetCollectionStartEnd(int g, int collectionSize,
                                              int &start, int &end) const {
  // default
  start = g * collectionSize / mGpuIds.size();
  end = (g+1) * collectionSize / mGpuIds.size();

  if (mGpuIds.size() == 1) {
    start = 0;
    end = collectionSize;
  }

  else if (mGpuIds.size() == 2) {
    if (collectionSize == 1) {
      if (g==0) {
        start = 0;
        end = 0;
      }
      else {
        start = 0;
        end = 1;
      }
    }

    else if (collectionSize == 9) {
      if (g==0) {
        start = 0;
        end = 4;
      }
      else {
        start = 4;
        end = 9;
      }
    }
  }
}
