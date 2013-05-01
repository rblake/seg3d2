
#include <map>
#include "sart.h"

using std::map;

void SARTContext::SetCurrentVolume(FloatVolumeType::Pointer vol) {
  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    mCurrentVolumeReconstruction[vi] = vol->GetBufferPointer()[vi];
  }
}

void SARTContext::SetCurrentVolume(float v) {
  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    mCurrentVolumeReconstruction[vi] = v; // 0.0384443121 * 100*100*100;
    
    int x,y,z;
    mGeometry.VolumeIndexToCoord(vi, x, y, z);
    if (z < 5)
      mCurrentVolumeReconstruction[vi] = 2.7e6;
  }
}

void SARTContext::SetBaselineProjection(FloatVolumeType::Pointer vol) {
  for (size_t pi=0; pi<mBaselineProjection.size(); pi++) {
    mBaselineProjection[pi] = vol->GetBufferPointer()[pi];
  }
}

void SARTContext::GetCurrentProjection(FloatVolumeType::Pointer vol) {
  FloatVolumeType::SizeType size;
  size.SetElement(0, mGeometry.GetDetectorSamplesWidth());
  size.SetElement(1, mGeometry.GetDetectorSamplesHeight());
  size.SetElement(2, mGeometry.GetNumProjectionAngles());
  vol->SetRegions(FloatVolumeType::RegionType(size));
  vol->Allocate();

  for (size_t pi=0; pi<mCurrentForwardProjection.size(); pi++) {
    vol->GetBufferPointer()[pi] = mCurrentForwardProjection[pi];
  }
}


void SARTContext::GetCurrentDeltaR(FloatVolumeType::Pointer vol) {
  FloatVolumeType::SizeType size;
  size.SetElement(0, mGeometry.GetDetectorSamplesWidth());
  size.SetElement(1, mGeometry.GetDetectorSamplesHeight());
  size.SetElement(2, mGeometry.GetNumProjectionAngles());
  vol->SetRegions(FloatVolumeType::RegionType(size));
  vol->Allocate();

  for (size_t pi=0; pi<mCurrentForwardProjection.size(); pi++) {
    vol->GetBufferPointer()[pi] = mCurrentForwardProjectionDeltaR[pi];
  }
}


void SARTContext::GetCurrentVolume(FloatVolumeType::Pointer vol) {
  FloatVolumeType::SizeType size;
  size.SetElement(0, mGeometry.GetVolumeSamplesX());
  size.SetElement(1, mGeometry.GetVolumeSamplesY());
  size.SetElement(2, mGeometry.GetVolumeSamplesZ());
  vol->SetRegions(FloatVolumeType::RegionType(size));
  vol->Allocate();

  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    vol->GetBufferPointer()[vi] = mCurrentVolumeReconstruction[vi];
  }
}


void SARTContext::CreateWeightsParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                        float stepSize, int samplesPerPixel) {

  vector<float> sourceSumSamples(mCurrentVolumeReconstruction.size());

  for (int p=id; p<mGeometry.GetNumProjectionAngles(); p+=nt) {
    std::cerr<<"create weights: angle "<<(p+1)<<" of "<<mGeometry.GetNumProjectionAngles()<<std::endl;

    const Vec3f sourcePosition = mGeometry.GetSourcePosition();


    // first integrate rays from the source through the volume
    for (int i=0; i<(int)mCurrentVolumeReconstruction.size(); i++) {
      sourceSumSamples[i] = 0;

      int volumeCoord[3];
      mGeometry.VolumeIndexToCoord(i, volumeCoord[0], volumeCoord[1], volumeCoord[2]);

      Vec3f fVolumeCoord(volumeCoord[0], volumeCoord[1], volumeCoord[2]);
      Vec3f worldCoord;
      mGeometry.VolumeToWorld(fVolumeCoord, worldCoord);

      // intersect ray starting at source position with volume bounding box
      Vec3f sourceRayDir = worldCoord - sourcePosition;
      float sourceDist = sourceRayDir.Length();
      sourceRayDir.Normalize();
      float smin, smax;
      if (!mGeometry.RayVolumeIntersection(sourcePosition, sourceRayDir, smin, smax))
        continue;
      smax = sourceDist;

      // integrate along the ray
      int nSourceSteps = (int)((smax-smin) / stepSize) + 1;
      float actualSourceStepSize = (smax-smin) / nSourceSteps;

      float sumSamples = 0;
      for (int ss=0; ss<nSourceSteps; ss++) {
        Vec3f sourceSamplePos = sourcePosition + ((ss+0.5f)*actualSourceStepSize + smin) * sourceRayDir;

        // trilinear sample at this position
        int volSourceIndices[8];
        float volSourceWeights[8];
        if (!mGeometry.WorldToTrilinearWeights(sourceSamplePos, volSourceIndices, volSourceWeights))
          continue;

        float valSource = 0;
        for (int i=0; i<8; i++)
          valSource += volSourceWeights[i] * mCurrentVolumeReconstruction[volSourceIndices[i]];

        sumSamples += actualSourceStepSize * valSource;
      }

      sourceSumSamples[i] = sumSamples;
    }


    for (int y=0; y<mGeometry.GetDetectorSamplesHeight(); y++) {
      for (int x=0; x<mGeometry.GetDetectorSamplesWidth(); x++) {
        int pi = mGeometry.ProjectionCoordToIndex(p, x, y);

        // map from voxels to weights for this detector pixel
        map<int,float> weightMap;


        for (int s=0; s<samplesPerPixel; s++) {

          // cast a ray out of this projection pixel
          Vec3f detectorRayOrigin, detectorRayDir;
          mGeometry.GetDetectorRay(p, x, y, (samplesPerPixel==1), detectorRayOrigin, detectorRayDir);

          // skip this detector pixel if no intersection
          float dmin, dmax;
          if (!mGeometry.RayVolumeIntersection(detectorRayOrigin, detectorRayDir, dmin, dmax))
            continue;
          if (dmin<0 || dmax<0)
            continue;

          // march along the ray from the detector
          int nDetectorSteps = (int)((dmax-dmin) / stepSize) + 1;
          float actualDetectorStepSize = (dmax-dmin) / nDetectorSteps;

          float sumDetectorSamples = 0;
          for (int ds=0; ds<nDetectorSteps; ds++) {
            Vec3f samplePos = detectorRayOrigin + ((ds+0.5f)*actualDetectorStepSize + dmin) * detectorRayDir;

            // trilinear sample at this position
            int volDetectorIndices[8];
            float volDetectorWeights[8];
            if (!mGeometry.WorldToTrilinearWeights(samplePos, volDetectorIndices, volDetectorWeights))
              continue;

            // source ray angular attenuation
            float sourceFrac = mGeometry.GetSourceIntensityThroughPoint(samplePos);
            if (sourceFrac == 0)
              continue;

            float valDetector = 0;
            for (int i=0; i<8; i++)
              valDetector += volDetectorWeights[i] * mCurrentVolumeReconstruction[volDetectorIndices[i]];

            sumDetectorSamples += actualDetectorStepSize * valDetector;

            // use precomputed source integrations
            Vec3f sourceRayDir = samplePos - sourcePosition;
            float sourceDist = sourceRayDir.Length();
            sourceRayDir.Normalize();

            float sumSourceSamples=0;
            for (int i=0; i<8; i++)
              sumSourceSamples += volDetectorWeights[i] * sourceSumSamples[volDetectorIndices[i]];

            // create a set of weights for this sample position
            float scatterAngle = (float)M_PI - acos(detectorRayDir.Dot(sourceRayDir));
            float sampleWeight = (actualDetectorStepSize *  // weight over the step length
                                  mFoamMaterial.GetScatterFactor(cos(scatterAngle)) * // scattering attenuation
                                  mFoamMaterial.GetMassAttenuationCoefficient() * // attenuation at the sample point
                                  exp(-mFoamMaterial.GetMassAttenuationCoefficient()*sumDetectorSamples) * // absorbption for detector ray
                                  exp(-mFoamMaterial.GetMassAttenuationCoefficient()*sumSourceSamples) * // absorbption for source ray
                                  sourceFrac * // attenuation from source inhomogeniety
                                  (1/(sourceDist * sourceDist))); // r^2 distance attenuation for point source

            // distribute this sample weight to the voxels around it
            for (int i=0; i<8; i++) {
              map<int,float>::iterator cur = weightMap.find(volDetectorIndices[i]);
              if (cur == weightMap.end())
                weightMap[volDetectorIndices[i]] = sampleWeight * volDetectorWeights[i] / samplesPerPixel;
              else
                cur->second += sampleWeight * volDetectorWeights[i] / samplesPerPixel;
            }
          }
        }


        // store all the weights from this detector ray
        ParallelEnterCriticalSection(mutex);
        for (map<int,float>::iterator wit=weightMap.begin(); wit!=weightMap.end(); ++wit) {
          int vi = wit->first;
          float w = wit->second;

          mProjectionToVolumeWeights[pi].push_back(pair<int,float>(vi,w));
          mVolumeToProjectionWeights[vi].push_back(pair<int,float>(pi,w));
        }
        ParallelLeaveCriticalSection(mutex);

      }
    }
  }

}


void SARTContext::CreateWeights(float stepSize, int samplesPerPixel, int numThreads) {

  mProjectionToVolumeWeights.clear();
  mProjectionToVolumeWeights.resize(mGeometry.GetTotalProjectionSamples());

  mVolumeToProjectionWeights.clear();
  mVolumeToProjectionWeights.resize(mGeometry.GetTotalVolumeSamples());

  ParallelExecutor(numThreads,
                   makeClassFunctor(this, &SARTContext::CreateWeightsParallel),
                   stepSize, samplesPerPixel);


  int totalWeights = 0;
  for (size_t pi=0; pi<mCurrentForwardProjection.size(); pi++) {
    totalWeights += mProjectionToVolumeWeights[pi].size();
  }
  printf("total weights created: %d\n", totalWeights);

  totalWeights = 0;
  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    totalWeights += mVolumeToProjectionWeights[vi].size();
  }
  printf("total weights created: %d\n", totalWeights);

}


void SARTContext::ForwardProject() {
  std::cerr<<"forward projection"<<std::endl;

  mCurrentForwardProjectionDeltaR.resize(mCurrentForwardProjection.size());

  for (size_t pi=0; pi<mCurrentForwardProjection.size(); pi++) {
    const vector< pair<int,float> > &weights = mProjectionToVolumeWeights[pi];

    double pv = 0;
    double pw = 0;

    for (size_t wi=0; wi<weights.size(); wi++) {
      const pair<int,float> &w = weights[wi];

      pv += (double)mCurrentVolumeReconstruction[w.first] * (double)w.second;
      pw += w.second;
    }

    if (pw>0) {
      mCurrentForwardProjection[pi] = pv;
      mCurrentForwardProjectionDeltaR[pi] = (mBaselineProjection[pi]-pv) / (pw);
    }
    else {
      mCurrentForwardProjection[pi] = 0;
      mCurrentForwardProjectionDeltaR[pi] = 0;
    }
  }
}


void SARTContext::BackwardProject(float lambda, int samplesPerPixel) {
  std::cerr<<"backward projection"<<std::endl;

  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    const vector< pair<int,float> > &weights = mVolumeToProjectionWeights[vi];

    double pv = 0;
    double pw = 0;

    for (size_t wi=0; wi<weights.size(); wi++) {
      const pair<int,float> &w = weights[wi];
      pv += (double)mCurrentForwardProjectionDeltaR[w.first] * (double)w.second;
      pw += w.second;
    }

    if (pw>0) {
      mCurrentVolumeReconstruction[vi] += lambda * pv/(pw);
      if (mCurrentVolumeReconstruction[vi] < 0)
        mCurrentVolumeReconstruction[vi] = 0;
    }
  }
}
