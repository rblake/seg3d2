
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

#include "markov.h"

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>

#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_functions.h>
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/helper_cuda.h>
#include <Application/BackscatterReconstruction/Algorithm/cuda_common/cutil_math.h>

// material volumes in the collection
float4 *cudaVolumeLinearCurrent = NULL;
float4 *cudaVolumeLinearCollection = NULL;
cudaArray *cudaVolumeArray00 = NULL;
cudaArray *cudaVolumeArray01 = NULL;
cudaArray *cudaVolumeArray02 = NULL;
cudaArray *cudaVolumeArray10 = NULL;
cudaArray *cudaVolumeArray11 = NULL;
cudaArray *cudaVolumeArray12 = NULL;
cudaArray *cudaVolumeArray20 = NULL;
cudaArray *cudaVolumeArray21 = NULL;
cudaArray *cudaVolumeArray22 = NULL;
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
cudaArray *cudaSourceArray00 = NULL;
cudaArray *cudaSourceArray01 = NULL;
cudaArray *cudaSourceArray02 = NULL;
cudaArray *cudaSourceArray10 = NULL;
cudaArray *cudaSourceArray11 = NULL;
cudaArray *cudaSourceArray12 = NULL;
cudaArray *cudaSourceArray20 = NULL;
cudaArray *cudaSourceArray21 = NULL;
cudaArray *cudaSourceArray22 = NULL;
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

// source ray info
float3 *cudaSourceRayVolOrigin = NULL;
float3 *cudaSourceRayVolDir = NULL;
float *cudaSourceRayTMin = NULL;
float *cudaSourceRayTMax = NULL;
float *cudaSourceScale = NULL; // accounts for source fall off and 1/r^2 attenuation
float *cudaSourceAttenuationCollection = NULL; // output

// detector ray info
float3 *cudaDetectorRayWorldOrigin = NULL;
float3 *cudaDetectorRayWorldDir = NULL;
float3 *cudaDetectorRayVolOrigin = NULL;
float3 *cudaDetectorRayVolDir = NULL;
float *cudaDetectorRayTMin = NULL;
float *cudaDetectorRayTMax = NULL;
float *cudaForwardProjectionCollection = NULL; // output

// info for calculating projection errors
float *cudaBaselineProjection = NULL;
float *cudaForwardProjectionError = NULL;


// info for optimizing which rays to cast
unsigned int *cudaRayIds = NULL;
float *cudaRayPriority = NULL;
float *cudaCurrentForwardProjection = NULL;


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
                               int collectionSize,
                               float3 matAtten,
                               float3 matDensity,
                               float voxelStepSize,
                               float *cudaSourceAttenuationCollection) {

  unsigned int gx = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int gy = blockIdx.y*blockDim.y + threadIdx.y;
  unsigned int gz = blockIdx.z*blockDim.z + threadIdx.z;
  if (gx>=volDim.x || gy>=volDim.y || gz>=volDim.z*collectionSize)
    return;
  unsigned int rayIndex = (gz%volDim.z)*volDim.x*volDim.y + gy*volDim.x + gx;
  int combo = gz/volDim.z;


  float3 volOrigin = cudaSourceRayVolOrigin[rayIndex];
  float3 volDir = cudaSourceRayVolDir[rayIndex];
  volDir = volDir / length(volDir);
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

  castSourceRays<<<dimGrid, dimBlock>>>(cudaSourceScale,
                                        cudaSourceRayVolOrigin,
                                        cudaSourceRayVolDir,
                                        cudaSourceRayTMin,
                                        cudaSourceRayTMax,
                                        volDim,
                                        collectionSize,
                                        matAtten,
                                        matDensity,
                                        mVoxelStepSize,
                                        cudaSourceAttenuationCollection);


  // bind output to source atten textures
  for (int c=0; c<collectionSize; c++) {
    cudaArray **cudaArray = NULL;
    switch (c) {
    case 0:
      cudaArray = &cudaSourceArray00;
      break;
    case 1:
      cudaArray = &cudaSourceArray01;
      break;
    case 2:
      cudaArray = &cudaSourceArray02;
      break;
    case 3:
      cudaArray = &cudaSourceArray10;
      break;
    case 4:
      cudaArray = &cudaSourceArray11;
      break;
    case 5:
      cudaArray = &cudaSourceArray12;
      break;
    case 6:
      cudaArray = &cudaSourceArray20;
      break;
    case 7:
      cudaArray = &cudaSourceArray21;
      break;
    case 8:
      cudaArray = &cudaSourceArray22;
      break;
    }


    cudaExtent volumeSize = make_cudaExtent(mGeometry.GetVolumeNodeSamplesX(),
                                            mGeometry.GetVolumeNodeSamplesY(),
                                            mGeometry.GetVolumeNodeSamplesZ());
    // copy data to 3D array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr   = make_cudaPitchedPtr((void *)&cudaSourceAttenuationCollection[mGeometry.GetTotalVolumeNodeSamples() * c],
                                              volumeSize.width*sizeof(float), volumeSize.width, volumeSize.height);
    copyParams.dstArray = *cudaArray;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyDeviceToDevice;
    checkCudaErrors(cudaMemcpy3DAsync(&copyParams));
  }
}


void MarkovContext::CudaGetSourceAttenuation(int collectionSize,
                                             vector< vector<float> > &sourceAttenuationCollection) const {
  // copy results back 
  vector<float> outputData(collectionSize * mGeometry.GetTotalVolumeNodeSamples());
  checkCudaErrors(cudaMemcpy(&outputData[0],
                             cudaSourceAttenuationCollection,
                             collectionSize * mGeometry.GetTotalVolumeNodeSamples()*sizeof(float),
                             cudaMemcpyDeviceToHost));

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
                                int3 detectorDim,
                                int3 volDim,
                                int collectionSize,
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

  float volLength = length(volP1-volP0);
  int numSteps = volLength / voxelStepSize + 1;


  int3 lastMat;
  float3 lastConc;
  float lastInterfaceT = tmin;

  float forwardProjection = 0;
  float sumDetectorAttenuation = 0;

  float lastInterfaceSourceAtten = sampleSourceAttenuation(combo, volOrigin + tmin*volDir, volDim);


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
                                    int3 detectorDim,
                                    int3 volDim,
                                    int collectionSize,
                                    float3 matAtten,
                                    float3 matDensity,
                                    float3 sourcePosition,
                                    float voxelStepSize,
                                    float *cudaForwardProjectionCollection) {

  unsigned int gx = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int gy = blockIdx.y*blockDim.y + threadIdx.y;
  unsigned int gz = blockIdx.z*blockDim.z + threadIdx.z;
  if (gx>=detectorDim.x || gy>=detectorDim.y || gz>=detectorDim.z*collectionSize)
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
                  detectorDim,
                  volDim,
                  collectionSize,
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
                                            int3 detectorDim,
                                            int3 volDim,
                                            int collectionSize,
                                            float3 matAtten,
                                            float3 matDensity,
                                            float3 sourcePosition,
                                            float voxelStepSize,
                                            float *cudaForwardProjectionCollection) {

  unsigned int gx = blockIdx.x*blockDim.x + threadIdx.x;
  if (gx>=detectorDim.x*detectorDim.y*detectorDim.z*collectionSize)
    return;

  //unsigned int rayIndexIndex = gx % (detectorDim.x*detectorDim.y*detectorDim.z);
  //int combo = gx / (detectorDim.x*detectorDim.y*detectorDim.z);
  unsigned int rayIndexIndex = gx / collectionSize;
  int combo = gx % collectionSize;

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
                    detectorDim,
                    volDim,
                    collectionSize,
                    matAtten,
                    matDensity,
                    sourcePosition,
                    voxelStepSize,
                    cudaForwardProjectionCollection);
  }
}




void MarkovContext::CudaForwardProject(int collectionSize) const {

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

  castAllDetectorRays<<<dimGrid, dimBlock>>>(cudaDetectorRayWorldOrigin,
                                             cudaDetectorRayWorldDir,
                                             cudaDetectorRayVolOrigin,
                                             cudaDetectorRayVolDir,
                                             cudaDetectorRayTMin,
                                             cudaDetectorRayTMax,
                                             detectorDim,
                                             volDim,
                                             collectionSize,
                                             matAtten,
                                             matDensity,
                                             sourcePosition,
                                             mVoxelStepSize,
                                             cudaForwardProjectionCollection);
}



void MarkovContext::CudaUpdateForwardProjection(int collectionSize,
                                                const Cone &attenChangeCone) const {

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
  prioritizeDetectorRays<<<dimGridPrioritize, dimBlockPrioritize>>>(cudaDetectorRayWorldOrigin,
                                                                    cudaDetectorRayWorldDir,
                                                                    cudaDetectorRayTMin,
                                                                    cudaDetectorRayTMax,
                                                                    mGeometry.GetTotalProjectionSamples(),
                                                                    make_float3(attenChangeCone.mOrigin[0],
                                                                                attenChangeCone.mOrigin[1],
                                                                                attenChangeCone.mOrigin[2]),
                                                                    make_float3(attenChangeCone.mDir[0],
                                                                                attenChangeCone.mDir[1],
                                                                                attenChangeCone.mDir[2]),
                                                                    attenChangeCone.mCosTheta,
                                                                    attenChangeCone.mMinDist,
                                                                    cudaRayIds,
                                                                    cudaRayPriority);


  // sort rays by priority
  thrust::device_ptr<unsigned int> thrustIds = thrust::device_pointer_cast(cudaRayIds);
  thrust::device_ptr<float> thrustPriorities = thrust::device_pointer_cast(cudaRayPriority);
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

  dim3 dimBlock(32,1,1);
  dim3 dimGrid(1+(collectionSize * mGeometry.GetTotalProjectionSamples()-1) / dimBlockPrioritize.x, 1, 1);


  castPrioritizedDetectorRays<<<dimGrid, dimBlock>>>(cudaRayIds,
                                                     cudaRayPriority,
                                                     cudaCurrentForwardProjection,
                                                     cudaDetectorRayWorldOrigin,
                                                     cudaDetectorRayWorldDir,
                                                     cudaDetectorRayVolOrigin,
                                                     cudaDetectorRayVolDir,
                                                     cudaDetectorRayTMin,
                                                     cudaDetectorRayTMax,
                                                     detectorDim,
                                                     volDim,
                                                     collectionSize,
                                                     matAtten,
                                                     matDensity,
                                                     sourcePosition,
                                                     mVoxelStepSize,
                                                     cudaForwardProjectionCollection);
}


void MarkovContext::CudaGetForwardProjection(int collectionSize,
                                             vector< vector<float> > &forwardProjectionCollection) const {
  // copy results back
  vector<float> outputData(collectionSize * mGeometry.GetTotalProjectionSamples());
  checkCudaErrors(cudaMemcpy(&outputData[0],
                             cudaForwardProjectionCollection,
                             collectionSize * mGeometry.GetTotalProjectionSamples()*sizeof(float),
                             cudaMemcpyDeviceToHost));

  forwardProjectionCollection.resize(collectionSize);
  for (int c=0; c<collectionSize; c++) {
    forwardProjectionCollection[c].resize(mGeometry.GetTotalProjectionSamples());
    for (int i=0; i<mGeometry.GetTotalProjectionSamples(); i++)
      forwardProjectionCollection[c][i] = outputData[c*mGeometry.GetTotalProjectionSamples() + i];
  }
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

  // compute squared errors
  int totalProjectionSamples = mGeometry.GetTotalProjectionSamples();
  dim3 dimBlock(32, 1, 1);
  dim3 dimGrid(1+((totalProjectionSamples*collectionSize)-1) / dimBlock.x, 1, 1);

  projectionToError<<<dimGrid, dimBlock>>>(totalProjectionSamples,
                                           collectionSize,
                                           cudaBaselineProjection,
                                           cudaForwardProjectionCollection,
                                           cudaForwardProjectionError);

  // use thrust to sum the errors for each material combo
  thrust::device_ptr<float> dev_ptr = thrust::device_pointer_cast(cudaForwardProjectionError);
  errors.resize(collectionSize);
  for (int c=0; c<collectionSize; c++) {
    errors[c] = (thrust::reduce(dev_ptr+c*totalProjectionSamples,
                                dev_ptr+(c+1)*totalProjectionSamples)
                 * mGeometry.GetDetectorPixelArea());
  }
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


void MarkovContext::CudaAcceptNextConfig(int c) const {
  // set the accepted volume as current
  checkCudaErrors(cudaMemcpyAsync(cudaVolumeLinearCurrent,
                                  cudaVolumeLinearCollection + c*mGeometry.GetTotalVolumeNodeSamples(),
                                  sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples(),
                                  cudaMemcpyDeviceToDevice));

  // set the accepted forward projection as current
  checkCudaErrors(cudaMemcpyAsync(cudaCurrentForwardProjection,
                                  cudaForwardProjectionCollection + c*mGeometry.GetTotalProjectionSamples(),
                                  sizeof(float)*mGeometry.GetTotalProjectionSamples(),
                                  cudaMemcpyDeviceToDevice));
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

  checkCudaErrors(cudaMemcpy(cudaVolumeLinearCurrent,
                             &volumeData[0],
                             sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples(),
                             cudaMemcpyHostToDevice));
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

  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
  cudaExtent volumeSize = make_cudaExtent(mGeometry.GetVolumeNodeSamplesX(),
                                          mGeometry.GetVolumeNodeSamplesY(),
                                          mGeometry.GetVolumeNodeSamplesZ());

  for (int c=0; c<NUM_MATERIALS*NUM_MATERIALS; c++) {

    float4 *cudaVolume = cudaVolumeLinearCollection + c*mGeometry.GetTotalVolumeNodeSamples();

    // copy current volume to collection volume
    checkCudaErrors(cudaMemcpyAsync(cudaVolume,
                                    cudaVolumeLinearCurrent,
                                    sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples(),
                                    cudaMemcpyDeviceToDevice));
  }

  for (int c=0; c<NUM_MATERIALS*NUM_MATERIALS; c++) {
    float4 *cudaVolume = cudaVolumeLinearCollection + c*mGeometry.GetTotalVolumeNodeSamples();

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

    
  for (int c=0; c<NUM_MATERIALS*NUM_MATERIALS; c++) {
    float4 *cudaVolume = cudaVolumeLinearCollection + c*mGeometry.GetTotalVolumeNodeSamples();

    cudaArray **cudaArray = NULL;
    switch (c) {
    case 0:
      cudaArray = &cudaVolumeArray00;
      break;
    case 1:
      cudaArray = &cudaVolumeArray01;
      break;
    case 2:
      cudaArray = &cudaVolumeArray02;
      break;
    case 3:
      cudaArray = &cudaVolumeArray10;
      break;
    case 4:
      cudaArray = &cudaVolumeArray11;
      break;
    case 5:
      cudaArray = &cudaVolumeArray12;
      break;
    case 6:
      cudaArray = &cudaVolumeArray20;
      break;
    case 7:
      cudaArray = &cudaVolumeArray21;
      break;
    case 8:
      cudaArray = &cudaVolumeArray22;
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



//==================================================================================================
//==================================================================================================
//==================================================================================================
void MarkovContext::CudaSetBaselineProjection() const {
  checkCudaErrors(cudaMalloc(&cudaBaselineProjection, sizeof(float)*mBaselineProjection.size()));
  checkCudaErrors(cudaMemcpy(cudaBaselineProjection, &mBaselineProjection[0], sizeof(float)*mBaselineProjection.size(), cudaMemcpyHostToDevice));
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

  // This will pick the best possible CUDA capable device
  int devID = findCudaDevice(0, NULL);

  // get device name
  cudaDeviceProp deviceProps;
  checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
  printf("CUDA device [%s]\n", deviceProps.name);


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
  checkCudaErrors(cudaMalloc(&cudaVolumeLinearCurrent, sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples()));
  checkCudaErrors(cudaMalloc(&cudaVolumeLinearCollection, sizeof(float4)*mGeometry.GetTotalVolumeNodeSamples()*NUM_MATERIALS*NUM_MATERIALS));
  cudaChannelFormatDesc channelDesc4 = cudaCreateChannelDesc<float4>();
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray00, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray01, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray02, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray10, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray11, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray12, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray20, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray21, &channelDesc4, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaVolumeArray22, &channelDesc4, volumeSize));

  SetTextureParams(&cudaVolumeTextures00);
  SetTextureParams(&cudaVolumeTextures01);
  SetTextureParams(&cudaVolumeTextures02);
  SetTextureParams(&cudaVolumeTextures10);
  SetTextureParams(&cudaVolumeTextures11);
  SetTextureParams(&cudaVolumeTextures12);
  SetTextureParams(&cudaVolumeTextures20);
  SetTextureParams(&cudaVolumeTextures21);
  SetTextureParams(&cudaVolumeTextures22);
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures00, cudaVolumeArray00, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures01, cudaVolumeArray01, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures02, cudaVolumeArray02, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures10, cudaVolumeArray10, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures11, cudaVolumeArray11, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures12, cudaVolumeArray12, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures20, cudaVolumeArray20, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures21, cudaVolumeArray21, channelDesc4));
  checkCudaErrors(cudaBindTextureToArray(cudaVolumeTextures22, cudaVolumeArray22, channelDesc4));


  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray00, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray01, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray02, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray10, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray11, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray12, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray20, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray21, &channelDesc, volumeSize));
  checkCudaErrors(cudaMalloc3DArray(&cudaSourceArray22, &channelDesc, volumeSize));

  SetTextureParams(&cudaSourceTextures00);
  SetTextureParams(&cudaSourceTextures01);
  SetTextureParams(&cudaSourceTextures02);
  SetTextureParams(&cudaSourceTextures10);
  SetTextureParams(&cudaSourceTextures11);
  SetTextureParams(&cudaSourceTextures12);
  SetTextureParams(&cudaSourceTextures20);
  SetTextureParams(&cudaSourceTextures21);
  SetTextureParams(&cudaSourceTextures22);
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures00, cudaSourceArray00, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures01, cudaSourceArray01, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures02, cudaSourceArray02, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures10, cudaSourceArray10, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures11, cudaSourceArray11, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures12, cudaSourceArray12, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures20, cudaSourceArray20, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures21, cudaSourceArray21, channelDesc));
  checkCudaErrors(cudaBindTextureToArray(cudaSourceTextures22, cudaSourceArray22, channelDesc));


  // upload all of the ray info
  checkCudaErrors(cudaMalloc(&cudaSourceRayVolOrigin, sizeof(float3)*mSourceRayVolOrigin.size()));
  checkCudaErrors(cudaMemcpy(cudaSourceRayVolOrigin, &mSourceRayVolOrigin[0], sizeof(float3)*mSourceRayVolOrigin.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaSourceRayVolDir, sizeof(float3)*mSourceRayVolDir.size()));
  checkCudaErrors(cudaMemcpy(cudaSourceRayVolDir, &mSourceRayVolDir[0], sizeof(float3)*mSourceRayVolDir.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaSourceRayTMin, sizeof(float)*mSourceRayTMin.size()));
  checkCudaErrors(cudaMemcpy(cudaSourceRayTMin, &mSourceRayTMin[0], sizeof(float)*mSourceRayTMin.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaSourceRayTMax, sizeof(float)*mSourceRayTMax.size()));
  checkCudaErrors(cudaMemcpy(cudaSourceRayTMax, &mSourceRayTMax[0], sizeof(float)*mSourceRayTMax.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaDetectorRayVolOrigin, sizeof(float3)*mDetectorRayVolOrigin.size()));
  checkCudaErrors(cudaMemcpy(cudaDetectorRayVolOrigin, &mDetectorRayVolOrigin[0], sizeof(float3)*mDetectorRayVolOrigin.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaDetectorRayVolDir, sizeof(float3)*mDetectorRayVolDir.size()));
  checkCudaErrors(cudaMemcpy(cudaDetectorRayVolDir, &mDetectorRayVolDir[0], sizeof(float3)*mDetectorRayVolDir.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaDetectorRayWorldOrigin, sizeof(float3)*mDetectorRayOrigin.size()));
  checkCudaErrors(cudaMemcpy(cudaDetectorRayWorldOrigin, &mDetectorRayOrigin[0], sizeof(float3)*mDetectorRayOrigin.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaDetectorRayWorldDir, sizeof(float3)*mDetectorRayDir.size()));
  checkCudaErrors(cudaMemcpy(cudaDetectorRayWorldDir, &mDetectorRayDir[0], sizeof(float3)*mDetectorRayDir.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaDetectorRayTMin, sizeof(float)*mDetectorRayTMin.size()));
  checkCudaErrors(cudaMemcpy(cudaDetectorRayTMin, &mDetectorRayTMin[0], sizeof(float)*mDetectorRayTMin.size(), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMalloc(&cudaDetectorRayTMax, sizeof(float)*mDetectorRayTMax.size()));
  checkCudaErrors(cudaMemcpy(cudaDetectorRayTMax, &mDetectorRayTMax[0], sizeof(float)*mDetectorRayTMax.size(), cudaMemcpyHostToDevice));

  checkCudaErrors(cudaMalloc(&cudaRayIds, sizeof(unsigned int)*mGeometry.GetTotalProjectionSamples()));
  checkCudaErrors(cudaMalloc(&cudaRayPriority, sizeof(float)*mGeometry.GetTotalProjectionSamples()));
  checkCudaErrors(cudaMalloc(&cudaCurrentForwardProjection, sizeof(float)*mGeometry.GetTotalProjectionSamples()));


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
    
    sourceScale[nvi] = mGeometry.GetSourceIntensityThroughPoint(voxelPosition) / (maxt*maxt);
  }

  checkCudaErrors(cudaMalloc(&cudaSourceScale, sizeof(float)*sourceScale.size()));
  checkCudaErrors(cudaMemcpy(cudaSourceScale, &sourceScale[0], sizeof(float)*sourceScale.size(), cudaMemcpyHostToDevice));


  checkCudaErrors(cudaMalloc(&cudaSourceAttenuationCollection, sizeof(float)*mGeometry.GetTotalVolumeNodeSamples() * NUM_MATERIALS*NUM_MATERIALS));
  checkCudaErrors(cudaMalloc(&cudaForwardProjectionCollection, sizeof(float)*mGeometry.GetTotalProjectionSamples() * NUM_MATERIALS*NUM_MATERIALS));
  checkCudaErrors(cudaMalloc(&cudaForwardProjectionError, sizeof(float)*mGeometry.GetTotalProjectionSamples() * NUM_MATERIALS*NUM_MATERIALS));

}


void MarkovContext::CudaShutdown() const {
  cudaDeviceReset();
}
