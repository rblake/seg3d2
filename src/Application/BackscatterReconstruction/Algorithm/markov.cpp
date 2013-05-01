
#include "markov.h"

#define USE_CUDA

MarkovContext::MarkovContext(const Geometry &g, const vector<Material> &materials,
                             int samplesPerPixel, float voxelStepSize, float energyRegularizationWeight)
  : mGeometry(g), 
    mMaterials(materials),
#ifdef USE_CUDA // single ray per pixel with cuda
    mSamplesPerPixel(1),
#else
    mSamplesPerPixel(samplesPerPixel),
#endif
    mVoxelStepSize(voxelStepSize),
    mEnergyRegularizationWeight(energyRegularizationWeight)
{

  mBaselineProjection.resize(mGeometry.GetTotalProjectionSamples(), 0);
  mCurrentForwardProjection.resize(mGeometry.GetTotalProjectionSamples(), 0);
  mCurrentVolumeReconstruction.resize(mGeometry.GetTotalVolumeNodeSamples(), 0);
  mCurrentVolumeSourceAttenuation.resize(mGeometry.GetTotalVolumeNodeSamples(), 0);


  // generate all detector rays that we'll use
  mDetectorRayOrigin.resize(mGeometry.GetTotalProjectionSamples() * samplesPerPixel);
  mDetectorRayDir.resize(mGeometry.GetTotalProjectionSamples() * samplesPerPixel);
  mDetectorRayVolOrigin.resize(mGeometry.GetTotalProjectionSamples() * samplesPerPixel);
  mDetectorRayVolDir.resize(mGeometry.GetTotalProjectionSamples() * samplesPerPixel);
  mDetectorRayTMin.resize(mGeometry.GetTotalProjectionSamples() * samplesPerPixel);
  mDetectorRayTMax.resize(mGeometry.GetTotalProjectionSamples() * samplesPerPixel);
  for (int pi=0; pi<mGeometry.GetTotalProjectionSamples(); pi++) {
    int p,x,y;
    mGeometry.ProjectionIndexToCoord(pi, p, x, y);
    
    for (int s=0; s<samplesPerPixel; s++) {
      int rayIndex = mGeometry.GetTotalProjectionSamples()*s + pi;
      mGeometry.GetDetectorRay(p, x, y, (mSamplesPerPixel==1), 
                               mDetectorRayOrigin[rayIndex],
                               mDetectorRayDir[rayIndex]);

      InitRay(mDetectorRayOrigin[rayIndex],
              mDetectorRayDir[rayIndex],
              mDetectorRayVolOrigin[rayIndex],
              mDetectorRayVolDir[rayIndex],
              mDetectorRayTMin[rayIndex],
              mDetectorRayTMax[rayIndex]);
    }
  }


  // generate all source rays that we'll use
  mSourceRayOrigin.resize(mGeometry.GetTotalVolumeNodeSamples());
  mSourceRayDir.resize(mGeometry.GetTotalVolumeNodeSamples());
  mSourceRayVolOrigin.resize(mGeometry.GetTotalVolumeNodeSamples());
  mSourceRayVolDir.resize(mGeometry.GetTotalVolumeNodeSamples());
  mSourceRayTMin.resize(mGeometry.GetTotalVolumeNodeSamples());
  mSourceRayTMax.resize(mGeometry.GetTotalVolumeNodeSamples());
  for (int nvi=0; nvi<mGeometry.GetTotalVolumeNodeSamples(); nvi++) {
    int x,y,z;
    mGeometry.VolumeIndexToNodeCoord(nvi, x,y,z);

    Vec3f voxelPosition;
    mGeometry.VolumeToWorld(Vec3f(x,y,z), voxelPosition);

    Vec3f diff = voxelPosition - mGeometry.GetSourcePosition();
    float tmax = diff.Length();
    Vec3f dir = diff / tmax;

    mSourceRayOrigin[nvi] = mGeometry.GetSourcePosition();
    mSourceRayDir[nvi] = dir;

    InitRay(mSourceRayOrigin[nvi],
            mSourceRayDir[nvi],
            mSourceRayVolOrigin[nvi],
            mSourceRayVolDir[nvi],
            mSourceRayTMin[nvi],
            mSourceRayTMax[nvi]);

    // source rays end at the node
    mSourceRayTMax[nvi] = std::min(mSourceRayTMax[nvi], tmax);
  }


#ifdef USE_CUDA
  CudaInitialize();
#endif
}


void MarkovContext::InitRay(const Vec3f &origin, const Vec3f &dir, 
                            Vec3f &volOrigin, Vec3f &volDir,
                            float &tmin, float &tmax) const {

  if (!mGeometry.RayVolumeIntersection(origin, dir, tmin, tmax)) {
    tmin = 0;
    tmax = -1;
    return;
  }

  mGeometry.WorldToVolume(origin, volOrigin);
  mGeometry.WorldToVolumeDir(dir, volDir);

  // make sure that when we round cast the volume ray position to int's, they are actually in the volume
  Vec3f volP0 = volOrigin + tmin*volDir;
  Vec3f volP1 = volOrigin + tmax*volDir;
  
  if (volP0[0] < 0 || volP0[1] < 0 || volP0[2] < 0 ||
      volP0[0] >= mGeometry.GetVolumeSamplesX() ||
      volP0[1] >= mGeometry.GetVolumeSamplesY() ||
      volP0[2] >= mGeometry.GetVolumeSamplesZ()) {
    tmin += 0.00001f;
  }

  if (volP1[0] < 0 || volP1[1] < 0 || volP1[2] < 0 ||
      volP1[0] >= mGeometry.GetVolumeSamplesX() ||
      volP1[1] >= mGeometry.GetVolumeSamplesY() ||
      volP1[2] >= mGeometry.GetVolumeSamplesZ()) {
    tmax -= 0.00001f;
  }

  volP0 = volOrigin + tmin*volDir;
  volP1 = volOrigin + tmax*volDir;
  if (volP0[0] < 0 || volP0[1] < 0 || volP0[2] < 0 ||
      volP0[0] >= mGeometry.GetVolumeSamplesX() ||
      volP0[1] >= mGeometry.GetVolumeSamplesY() ||
      volP0[2] >= mGeometry.GetVolumeSamplesZ() ||
      volP1[0] < 0 || volP1[1] < 0 || volP1[2] < 0 ||
      volP1[0] >= mGeometry.GetVolumeSamplesX() ||
      volP1[1] >= mGeometry.GetVolumeSamplesY() ||
      volP1[2] >= mGeometry.GetVolumeSamplesZ()) {
    tmin = 0;
    tmax = -1;
  }

}


void MarkovContext::SetCurrentVolume(ByteVolumeType::Pointer vol) {
  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    mCurrentVolumeReconstruction[vi] = vol->GetBufferPointer()[vi];
  }
}

void MarkovContext::SetCurrentVolume(unsigned char v) {
  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    mCurrentVolumeReconstruction[vi] = v;

    int x,y,z;
    mGeometry.VolumeIndexToNodeCoord(vi, x, y, z);
    if (z < 1)
      mCurrentVolumeReconstruction[vi] = 2;

    //if (z > 4)
    //mCurrentVolumeReconstruction[vi] = 0;
    
  }
}

void MarkovContext::GetCurrentVolume(ByteVolumeType::Pointer vol) {
  ByteVolumeType::SizeType size;
  size.SetElement(0, mGeometry.GetVolumeNodeSamplesX());
  size.SetElement(1, mGeometry.GetVolumeNodeSamplesY());
  size.SetElement(2, mGeometry.GetVolumeNodeSamplesZ());
  vol->SetRegions(ByteVolumeType::RegionType(size));
  vol->Allocate();

  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {
    vol->GetBufferPointer()[vi] = mCurrentVolumeReconstruction[vi];
  }
}



void MarkovContext::SetBaselineProjection(FloatVolumeType::Pointer vol) {
  for (size_t pi=0; pi<mBaselineProjection.size(); pi++) {
    int p,x,y;
    mGeometry.ProjectionIndexToCoord(pi, p, x, y);
    float val = vol->GetBufferPointer()[pi];

    // gain/offset correction
    val = (val - mGeometry.GetDetectorOffset()) / mGeometry.GetDetectorGain();

    mBaselineProjection[pi] = val;
  }

#ifdef USE_CUDA
  CudaSetBaselineProjection();
#endif
}


void MarkovContext::GetCurrentProjection(FloatVolumeType::Pointer vol) {
  FloatVolumeType::SizeType size;
  size.SetElement(0, mGeometry.GetDetectorSamplesWidth());
  size.SetElement(1, mGeometry.GetDetectorSamplesHeight());
  size.SetElement(2, mGeometry.GetNumProjectionAngles());
  vol->SetRegions(FloatVolumeType::RegionType(size));
  vol->Allocate();

  for (size_t pi=0; pi<mCurrentForwardProjection.size(); pi++) {
    vol->GetBufferPointer()[pi] = mCurrentForwardProjection[pi] * mGeometry.GetDetectorGain() + mGeometry.GetDetectorOffset();
  }
}


void MarkovContext::GetCurrentVolumeSourceAttenuation(FloatVolumeType::Pointer vol) const {
  FloatVolumeType::SizeType size;

  size.SetElement(0, mGeometry.GetVolumeNodeSamplesX());
  size.SetElement(1, mGeometry.GetVolumeNodeSamplesY());
  size.SetElement(2, mGeometry.GetVolumeNodeSamplesZ());

  vol->SetRegions(FloatVolumeType::RegionType(size));
  vol->Allocate();

  for (size_t vi=0; vi<mCurrentVolumeSourceAttenuation.size(); vi++) {
    vol->GetBufferPointer()[vi] = mCurrentVolumeSourceAttenuation[vi];
  }
}


void MarkovContext::FindInitialGuess() {

  float best_v = 1e10;
  int best_al = -1;
  int best_foam = -1;

  for (int al=0; al<mGeometry.GetVolumeNodeSamplesZ()-2; al++) {
    for (int foam=al+1; foam<mGeometry.GetVolumeNodeSamplesZ(); foam++) {

      for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {

        int x,y,z;
        mGeometry.VolumeIndexToNodeCoord(vi, x, y, z);

        if (z < al)
          mCurrentVolumeReconstruction[vi] = 0;

        else if (z == al)
          mCurrentVolumeReconstruction[vi] = 2;

        else if (z < foam)
          mCurrentVolumeReconstruction[vi] = 1;

        else
          mCurrentVolumeReconstruction[vi] = 0;
      }

      ComputeSourceAttenuation();
      ComputeForwardProjection();
      double err = ComputeTotalError(mCurrentVolumeReconstruction, mCurrentForwardProjection);

      std::cerr<<al<<" "<<foam<<" "<<err<<std::endl;

      if (err < best_v) {
        best_v = err;
        best_al = al;
        best_foam = foam;
      }

    }
  }

  std::cerr<<"best: "<<best_al<<" "<<best_foam<<" "<<best_v<<std::endl;

  // reset to best combination
  for (size_t vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {

    int x,y,z;
    mGeometry.VolumeIndexToNodeCoord(vi, x, y, z);

    if (z < best_al)
      mCurrentVolumeReconstruction[vi] = 0;

    else if (z == best_al)
      mCurrentVolumeReconstruction[vi] = 2;

    else if (z < best_foam)
      mCurrentVolumeReconstruction[vi] = 1;

    else
      mCurrentVolumeReconstruction[vi] = 0;
  }

}


void LerpBitMatCollection(int vals[NUM_MATERIALS*NUM_MATERIALS][2][2][2], Vec3f pos, int collectionSize,
                          float concentrations[NUM_MATERIALS*NUM_MATERIALS][NUM_MATERIALS]) {

  float factor[8] = {
    ((1-pos[0])*(1-pos[1])*(1-pos[2])),
    ((1-pos[0])*(1-pos[1])*(  pos[2])),
    ((1-pos[0])*(  pos[1])*(1-pos[2])),
    ((1-pos[0])*(  pos[1])*(  pos[2])),
    ((  pos[0])*(1-pos[1])*(1-pos[2])),
    ((  pos[0])*(1-pos[1])*(  pos[2])),
    ((  pos[0])*(  pos[1])*(1-pos[2])),
    ((  pos[0])*(  pos[1])*(  pos[2]))
  };


  for (int c=0; c<collectionSize; c++) {
    for (int m=0; m<NUM_MATERIALS; m++) {

      concentrations[c][m] = (((vals[c][0][0][0]>>m)&1) * factor[0] +
                              ((vals[c][0][0][1]>>m)&1) * factor[1] +
                              ((vals[c][0][1][0]>>m)&1) * factor[2] +
                              ((vals[c][0][1][1]>>m)&1) * factor[3] +
                              ((vals[c][1][0][0]>>m)&1) * factor[4] +
                              ((vals[c][1][0][1]>>m)&1) * factor[5] +
                              ((vals[c][1][1][0]>>m)&1) * factor[6] +
                              ((vals[c][1][1][1]>>m)&1) * factor[7]);
    }
  }
}



// get a list of intervals corresponding to marching the given ray through the volume
bool MarkovContext::GetRayVoxelIntersections(const vector< vector<unsigned char> > &volumeReconstructionCollection, 
                                             int collectionSize,
                                             int rayIndex, bool detectorRay,
                                             float tmin, float tmax,
                                             vector<float> &ts, vector< vector<unsigned char> > &materials) const {

  // clear output
  ts.clear();
  materials.resize(collectionSize);
  for (int i=0; i<materials.size(); i++)
    materials[i].clear();

  float _tmin = detectorRay ? mDetectorRayTMin[rayIndex] : mSourceRayTMin[rayIndex];
  float _tmax = detectorRay ? mDetectorRayTMax[rayIndex] : mSourceRayTMax[rayIndex];
  tmin = std::max(tmin, _tmin);
  tmax = std::min(tmax, _tmax);
  if (tmin > tmax)
    return false;


  const Vec3f &volOrigin = detectorRay ? mDetectorRayVolOrigin[rayIndex] : mSourceRayVolOrigin[rayIndex];
  const Vec3f &volDir = detectorRay ? mDetectorRayVolDir[rayIndex] : mSourceRayVolDir[rayIndex];

  Vec3f volP0 = volOrigin + tmin*volDir;
  Vec3f volP1 = volOrigin + tmax*volDir;


#if 0 // should never happen!
  if (volP0[0] < 0 || volP0[1] < 0 || volP0[2] < 0 ||
      volP0[0] >= mGeometry.GetVolumeSamplesX() ||
      volP0[1] >= mGeometry.GetVolumeSamplesY() ||
      volP0[2] >= mGeometry.GetVolumeSamplesZ() ||
      volP1[0] < 0 || volP1[1] < 0 || volP1[2] < 0 ||
      volP1[0] >= mGeometry.GetVolumeSamplesX() ||
      volP1[1] >= mGeometry.GetVolumeSamplesY() ||
      volP1[2] >= mGeometry.GetVolumeSamplesZ()) {
    std::cerr<<"ray outside of volume!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    return false; // ray was just grazing and adjust start and end to be inside volume moved it back out
  }
#endif


  float volLength = (volP1-volP0).Length();
  float worldLength = tmax-tmin;

  int numSteps = volLength / mVoxelStepSize + 1;
  float volStepLength = volLength / numSteps;
  float worldStepLength = worldLength / volLength;

  unsigned char lastMat[NUM_MATERIALS*NUM_MATERIALS];
  for (int i=0; i<NUM_MATERIALS*NUM_MATERIALS; i++)
    lastMat[i] = 0xFF;

  float lastMatConcentrations[NUM_MATERIALS*NUM_MATERIALS][NUM_MATERIALS];
  for (int i=0; i<NUM_MATERIALS*NUM_MATERIALS; i++) {
    for (int j=0; j<NUM_MATERIALS; j++) {
      lastMatConcentrations[i][j] = 0;
    }
  }


  ts.push_back(tmin);

  int cellMats[NUM_MATERIALS*NUM_MATERIALS][2][2][2];
  Vec3i lastVolPi(-1,-1,-1);


  for (int step=0; step<numSteps; step++) {
    float fstep = (step+0.5f) / numSteps;
    Vec3f volP = (1-fstep)*volP0 + fstep*volP1;

    Vec3i volPi((int)volP[0], (int)volP[1], (int)volP[2]);
    Vec3f volPf(volP[0]-volPi[0], volP[1]-volPi[1], volP[2]-volPi[2]);

    // only load new cell mats if we're in a different cell
    if (lastVolPi[0]!=volPi[0] || lastVolPi[1]!=volPi[1] || lastVolPi[2]!=volPi[2]) {
      for (int z=0; z<2; z++) {
        for (int y=0; y<2; y++) {
          for (int x=0; x<2; x++) {
            int vxyz = mGeometry.VolumeNodeCoordToIndex(volPi[0]+x, volPi[1]+y, volPi[2]+z);
            for (int c=0; c<collectionSize; c++) {
              cellMats[c][x][y][z] = 1<<volumeReconstructionCollection[c][vxyz];
            }
          }
        }
      }

      lastVolPi = volPi;
    }

    // get new interpolated material concentrations
    float thisMatConcentrations[NUM_MATERIALS*NUM_MATERIALS][NUM_MATERIALS];
    LerpBitMatCollection(cellMats, volPf, collectionSize, thisMatConcentrations);



    // find the current material for each combination
    unsigned char thisMat[NUM_MATERIALS*NUM_MATERIALS];
    for (int c=0; c<collectionSize; c++) {
      thisMat[c] = 0;
      float thisMatVal = thisMatConcentrations[c][0];

      for (int m=1; m<mMaterials.size(); m++) {
        float matVal = thisMatConcentrations[c][m];
        if (matVal > thisMatVal) {
          thisMatVal = matVal;
          thisMat[c] = m;
        }
      }
    }


    // mark a material change if any of the combinations passed an interface
    if (step>0) {

      // assume linear concentration changes - interpolate the step position
      float avet = 0;
      int nt=0;
      for (int c=0; c<collectionSize; c++) {
        if (thisMat[c]!=lastMat[c]) {
          float t1MatConc = thisMatConcentrations[c][thisMat[c]];
          float t1OtherConc = 0;
          for (int i=0; i<NUM_MATERIALS-1; i++) {
            t1OtherConc += thisMatConcentrations[c][(thisMat[c]+1+i)%NUM_MATERIALS];
          }

          float t0MatConc = lastMatConcentrations[c][lastMat[c]];
          float t0OtherConc = 0;
          for (int i=0; i<NUM_MATERIALS-1; i++) {
            t0OtherConc += lastMatConcentrations[c][(lastMat[c]+1+i)%NUM_MATERIALS];
          }

          float tx = (t0MatConc-t0OtherConc) / ((t1MatConc-t0OtherConc) - (t1OtherConc-t0MatConc));

          if (tx<0) tx=0;
          if (tx>1) tx=1;

          avet += tx;
          nt++;
        }
      }



      if (nt>0) {

        avet /= nt;
        float t0 = tmax * ((float)(step-0.5f)/numSteps) + tmin * (1-(float)(step-0.5f)/numSteps);
        float t1 = tmax * ((float)(step+0.5f)/numSteps) + tmin * (1-(float)(step+0.5f)/numSteps);
        float thist = (1-avet)*t0 + avet*t1;

        // add non-interface steps so we cap our max step size at ~1 voxel (frequency of source atten)
        float lastt = ts.back();
        float thisStepLen = thist - lastt;
        int substeps = (int)(thisStepLen / worldStepLength);
        for (int i=0; i<substeps; i++) {
          ts.push_back(lastt + ((i+0.5f)/substeps) * (thist-lastt));
          for (int c=0; c<collectionSize; c++) {
            materials[c].push_back(lastMat[c]);
          }
        }

        ts.push_back(thist);
        for (int c=0; c<collectionSize; c++) {
          materials[c].push_back(lastMat[c]);
        }
      }
      /*
      else {
        avet = 0.5f;
        float t0 = tmax * ((float)(step-0.5f)/numSteps) + tmin * (1-(float)(step-0.5f)/numSteps);
        float t1 = tmax * ((float)(step+0.5f)/numSteps) + tmin * (1-(float)(step+0.5f)/numSteps);
        float thist = (1-avet)*t0 + avet*t1;

        // add non-interface steps so we cap our max step size at ~1 voxel (frequency of source atten)
        float lastt = ts.back();
        float thisStepLen = thist - lastt;
        int substeps = (int)(thisStepLen / worldStepLength);
        for (int i=0; i<substeps; i++) {
          ts.push_back(lastt + ((i+0.5f)/substeps) * (thist-lastt));
          for (int c=0; c<collectionSize; c++) {
            materials[c].push_back(lastMat[c]);
          }
        }

        ts.push_back(thist);
        for (int c=0; c<collectionSize; c++) {
          materials[c].push_back(lastMat[c]);
        }
      }
      */
    }

    // update last material concentrations
    for (int c=0; c<collectionSize; c++) {
      for (int j=0; j<NUM_MATERIALS; j++) {
        lastMatConcentrations[c][j] = thisMatConcentrations[c][j];
      }
      lastMat[c] = thisMat[c];
    }
  }

  // add final material s
  float lastt = ts.back();
  float thisStepLen = tmax - lastt;
  int substeps = (int)(thisStepLen / worldStepLength);
  for (int i=0; i<substeps; i++) {
    ts.push_back(lastt + ((i+0.5f)/substeps) * (tmax-lastt));
    for (int c=0; c<collectionSize; c++) {
      materials[c].push_back(lastMat[c]);
    }
  }

  ts.push_back(tmax);
  for (int c=0; c<collectionSize; c++) {
    materials[c].push_back(lastMat[c]);
  }


  return true;
}


void MarkovContext::GetAttenChangeCone(const Cone &changeCone, Cone &attenChangeCone) const {

  attenChangeCone = changeCone;

  for (int nvi=0; nvi<mGeometry.GetTotalVolumeNodeSamples(); nvi++) {
    int x,y,z;
    mGeometry.VolumeIndexToNodeCoord(nvi, x,y,z);

    Vec3f voxelPosition;
    mGeometry.VolumeToWorld(Vec3f(x,y,z), voxelPosition);

    if (changeCone.ContainsPoint(voxelPosition)) {
      Vec3f voxelMin, voxelMax;
      mGeometry.VolumeToWorld(Vec3f(x-1,y-1,z-1), voxelMin);
      mGeometry.VolumeToWorld(Vec3f(x+1,y+1,z+1), voxelMax);
      attenChangeCone.Expand(voxelMin, voxelMax);
    }
  }
}


// integrate the attenuation coefficients from the source into the volume
void MarkovContext::UpdateSourceAttenuationCollectionSample(const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                                            const vector<float> &currentVolumeSourceAttenuation,
                                                            const Cone *sourceChangeCone,
                                                            int nvi,
                                                            Cone &attenChangeCone,
                                                            vector< vector<float> > &sourceAttenuation) const {

  Vec3f sourcePosition = mGeometry.GetSourcePosition();
  const int collectionSize = volumeReconstructionCollection.size();

  vector<float> ts;
  vector< vector<unsigned char> > mats;
  vector<float> matlens(mMaterials.size());

  int x,y,z;
  mGeometry.VolumeIndexToNodeCoord(nvi, x,y,z);

  Vec3f voxelPosition;
  mGeometry.VolumeToWorld(Vec3f(x,y,z), voxelPosition);


  // the ray to the voxel doesn't go through the change
  if (sourceChangeCone && !sourceChangeCone->ContainsPoint(voxelPosition)) {
    for (int c=0; c<collectionSize; c++)
      sourceAttenuation[c][nvi] = currentVolumeSourceAttenuation[nvi];
  }

  // the ray to the voxel goes through the change - full ray trace
  else {
    Vec3f voxelMin, voxelMax;
    mGeometry.VolumeToWorld(Vec3f(x-1,y-1,z-1), voxelMin);
    mGeometry.VolumeToWorld(Vec3f(x+1,y+1,z+1), voxelMax);
    attenChangeCone.Expand(voxelMin, voxelMax);

    for (int c=0; c<collectionSize; c++)
      sourceAttenuation[c][nvi] = 1;

    Vec3f diff = voxelPosition - sourcePosition;
    float maxt = diff.Length();
    Vec3f dir = diff / maxt;

    if (!GetRayVoxelIntersections(volumeReconstructionCollection, (int)volumeReconstructionCollection.size(),
                                  nvi, false,
                                  0, maxt, ts, mats)) {
      //std::cerr<<"no ray intersections for voxel???"<<std::endl;
    }
    else {

      // attenuate for each material
      for (int c=0; c<collectionSize; c++) {

        for (int m=0; m<matlens.size(); m++)
          matlens[m] = 0;

        const vector<unsigned char> &tmats = mats[c];
        for (int i=0; i<tmats.size(); i++) {
          matlens[tmats[i]] += ts[i+1]-ts[i];
        }

        double atten = 1;
        for (int i=0; i<(int)matlens.size(); i++) {
          atten *= exp(-matlens[i] * (mMaterials[i].GetDensity() * mMaterials[i].GetMassAttenuationCoefficient()));
        }
        sourceAttenuation[c][nvi] = (float)atten;
      }
    }

    // 1/r^2 attenuation, radial source falloft
    for (int c=0; c<collectionSize; c++) {
      sourceAttenuation[c][nvi] *= mGeometry.GetSourceIntensityThroughPoint(voxelPosition) / (maxt*maxt);
    }

  }
}


void MarkovContext::ComputeSourceAttenuation() {

  vector< vector<float> > sourceAttenuationCollection;
  sourceAttenuationCollection.push_back(mCurrentVolumeSourceAttenuation);


#ifndef USE_CUDA
  vector< vector<unsigned char> > volumeReconstructionCollection;
  volumeReconstructionCollection.push_back(mCurrentVolumeReconstruction);

  vector<float> bogusCurrentVolumeSourceAttenuation;
  Cone attenChangeCone(mGeometry.GetSourcePosition(), Vec3f(0,0,-1), 0);
  
  for (int nvi=0; nvi<(int)mCurrentVolumeSourceAttenuation.size(); nvi++) {
    UpdateSourceAttenuationCollectionSample(volumeReconstructionCollection,
                                            bogusCurrentVolumeSourceAttenuation,
                                            NULL,
                                            nvi,
                                            attenChangeCone,
                                            sourceAttenuationCollection);
  }

#else

  CudaSetCurrentVolume(mCurrentVolumeReconstruction);
  CudaSetVolumeCollection(GibbsProposal(-1,-1));

  // test that we set the volumes correctly
  vector< vector<unsigned char> > volCollection;
  CudaGetVolumeCollection(volCollection);
  for (int c=0; c<volCollection.size(); c++) {
    for (int i=0; i<volCollection[c].size(); i++) {
      if (volCollection[c][i] != mCurrentVolumeReconstruction[i]) {
        std::cerr<<"volume collection not set correctly!"<<std::endl;
      }
    }
  }


  CudaComputeSourceAttenuation(1, &sourceAttenuationCollection);

#endif


  mCurrentVolumeSourceAttenuation = sourceAttenuationCollection[0];

}


void MarkovContext::UpdateSourceAttenuationCollectionParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                                              const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                                              const vector<float> &currentVolumeSourceAttenuation,
                                                              const Cone &sourceChangeCone,
                                                              Cone &attenChangeCone,
                                                              vector< vector<float> > &sourceAttenuation) const {

  ParallelEnterCriticalSection(mutex);
  Cone myAttenChangeCone = attenChangeCone;
  ParallelLeaveCriticalSection(mutex);


  for (int nvi=id; nvi<(int)currentVolumeSourceAttenuation.size(); nvi+=nt) {
    UpdateSourceAttenuationCollectionSample(volumeReconstructionCollection,
                                            currentVolumeSourceAttenuation,
                                            &sourceChangeCone,
                                            nvi,
                                            myAttenChangeCone,
                                            sourceAttenuation);
  }

  ParallelEnterCriticalSection(mutex);
  attenChangeCone.MergeSameAxisOrigin(myAttenChangeCone);
  ParallelLeaveCriticalSection(mutex);
}


void MarkovContext::LerpNodeValuesCollection(const vector< vector<float> > &valueCollection,
                                             int collectionSize,
                                             const Vec3f &volumePos,
                                             vector<float> &lerpedValues) const {

  Vec3<int> samplePosI(std::max(0, std::min(mGeometry.GetVolumeSamplesX()-1, (int)volumePos[0])),
                       std::max(0, std::min(mGeometry.GetVolumeSamplesY()-1, (int)volumePos[1])),
                       std::max(0, std::min(mGeometry.GetVolumeSamplesZ()-1, (int)volumePos[2])));

  Vec3f samplePosF(volumePos[0] - samplePosI[0],
                   volumePos[1] - samplePosI[1],
                   volumePos[2] - samplePosI[2]);

  int indices[8] = {
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+0, samplePosI[1]+0, samplePosI[2]+0),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+0, samplePosI[1]+0, samplePosI[2]+1),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+0, samplePosI[1]+1, samplePosI[2]+0),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+0, samplePosI[1]+1, samplePosI[2]+1),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+1, samplePosI[1]+0, samplePosI[2]+0),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+1, samplePosI[1]+0, samplePosI[2]+1),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+1, samplePosI[1]+1, samplePosI[2]+0),
    mGeometry.VolumeNodeCoordToIndex(samplePosI[0]+1, samplePosI[1]+1, samplePosI[2]+1)
  };

  float factors[8] = {
    ((1-samplePosF[0])*(1-samplePosF[1])*(1-samplePosF[2])),
    ((1-samplePosF[0])*(1-samplePosF[1])*(  samplePosF[2])),
    ((1-samplePosF[0])*(  samplePosF[1])*(1-samplePosF[2])),
    ((1-samplePosF[0])*(  samplePosF[1])*(  samplePosF[2])),
    ((  samplePosF[0])*(1-samplePosF[1])*(1-samplePosF[2])),
    ((  samplePosF[0])*(1-samplePosF[1])*(  samplePosF[2])),
    ((  samplePosF[0])*(  samplePosF[1])*(1-samplePosF[2])),
    ((  samplePosF[0])*(  samplePosF[1])*(  samplePosF[2]))
  };

  lerpedValues.resize(valueCollection.size());
  for (int i=0; i<collectionSize; i++) {
    //const vector<float> &values = valueCollection[i];
    const float *values = &valueCollection[i][0];
    lerpedValues[i] = (values[indices[0]]*factors[0] + 
                       values[indices[1]]*factors[1] + 
                       values[indices[2]]*factors[2] + 
                       values[indices[3]]*factors[3] + 
                       values[indices[4]]*factors[4] + 
                       values[indices[5]]*factors[5] + 
                       values[indices[6]]*factors[6] + 
                       values[indices[7]]*factors[7]);
  }
}



void MarkovContext::GetForwardProjectionCollectionForPixel(int p, int x, int y,
                                                           int rayIndex,
                                                           float t0, float t1,
                                                           int collectionSize,
                                                           const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                                           const vector< vector<float> > &sourceAttenuationCollection,
                                                           vector<float> &forwardProjectionPixelCollection,
                                                           vector<float> &ts,
                                                           vector< vector<unsigned char> > &mats,
                                                           vector<float> startSourceAttenCollection,
                                                           vector<float> endSourceAttenCollection,
                                                           vector<float> &sumDetectorAttenuation) const {
  

  Vec3f sourcePosition = mGeometry.GetSourcePosition();

  forwardProjectionPixelCollection.resize(collectionSize);
  sumDetectorAttenuation.resize(collectionSize);
  for (int c=0; c<collectionSize; c++) {
    forwardProjectionPixelCollection[c] = 0;
    sumDetectorAttenuation[c] = 0;
  }


  if (!GetRayVoxelIntersections(volumeReconstructionCollection, collectionSize,
                                rayIndex, true,
                                t0, t1, ts, mats)) {
    //std::cerr<<"no ray intersections for voxel???"<<std::endl;
  }

  else {
    const Vec3f &detectorRayOrigin = mDetectorRayOrigin[rayIndex];
    const Vec3f &detectorRayDir = mDetectorRayDir[rayIndex];

    Vec3f startSamplePos = detectorRayOrigin + ts.front()*detectorRayDir;
    Vec3f startSamplePosV;
    mGeometry.WorldToVolume(startSamplePos, startSamplePosV);

    LerpNodeValuesCollection(sourceAttenuationCollection, collectionSize, startSamplePosV, startSourceAttenCollection);
    
    for (int v=0; v<ts.size()-1; v++) {

      float tdist = ts[v+1]-ts[v];
      float tcenter = 0.5 * (ts[v+1]+ts[v]);

      Vec3f samplePos = detectorRayOrigin + tcenter*detectorRayDir;
      Vec3f sourceRayDir = samplePos - sourcePosition;
      float sourceDist = sourceRayDir.Length();
      sourceRayDir /= sourceDist;
      float ncosScatterAngle = detectorRayDir.Dot(sourceRayDir);


      Vec3f endSamplePos = detectorRayOrigin + ts[v+1]*detectorRayDir;
      Vec3f endSamplePosV;
      mGeometry.WorldToVolume(endSamplePos, endSamplePosV);

      LerpNodeValuesCollection(sourceAttenuationCollection, collectionSize, endSamplePosV, endSourceAttenCollection);

      // all material configurations are different
      for (int c=0; c<collectionSize; c++) {
        int materiali = mats[c][v];
        const Material &material = mMaterials[materiali];

        float voxelAttenuation = -tdist * material.GetDensity() * material.GetMassAttenuationCoefficient();
        float endSourceAtten = endSourceAttenCollection[c];
        float &startSourceAtten = startSourceAttenCollection[c];

        // attenuationFactor = int_0^1 [ ((1-t)*startSourceAtten + t*endSourceAtten) * exp(t* (-tdist*density*atten)) ] dt
        // a = startSourceAtten
        // b = endSourceAtten
        // c = voxelAttenuation
        // integral = (e^c * (a+b(c-1)) - a(c+1)+b) / c^2
        float attenuationFactor;
        if (voxelAttenuation == 0) {
          attenuationFactor = (startSourceAtten+endSourceAtten)*0.5f;
        }
        else {
          double a = startSourceAtten;
          double b = endSourceAtten;
          double c = voxelAttenuation;
          attenuationFactor = (exp(c) * (a+b*(c-1)) - a*(c+1)+b) / (c*c);
        }

        // attenuation between detector and start
        attenuationFactor *= exp(sumDetectorAttenuation[c]);

        /*
        // monte carlo integration to double check it
        {
          float mcAttenuationFactor = 0;
          for (int i=0; i<1000; i++) {
            float f = (i+0.5f)/1000;
            mcAttenuationFactor += (startSourceAtten*(1-f) + endSourceAtten*f) * exp(f * -tdist * material.GetDensity() * material.GetMassAttenuationCoefficient());
          }
          mcAttenuationFactor /= 1000;
          mcAttenuationFactor *= exp(sumDetectorAttenuation[c]);

          if (attenuationFactor / mcAttenuationFactor > (1+1e-4) ||
              mcAttenuationFactor / attenuationFactor > (1+1e-4))
            {
              std::cerr<<attenuationFactor / mcAttenuationFactor<<std::endl;
              std::cerr<<attenuationFactor<<" "<<mcAttenuationFactor<<std::endl;
            }
        }
        */

        forwardProjectionPixelCollection[c] += (tdist *
                                                material.GetScatterFactor(0.5 * (1+ncosScatterAngle)) * // scattering attenuation
                                                attenuationFactor);

        // accumulate attenuation for detector ray
        sumDetectorAttenuation[c] += voxelAttenuation;

        startSourceAtten = endSourceAtten;
      }
    }
  }
}


// use the current volume guess and source attenuation to do a forward projection
void MarkovContext::ComputeForwardProjection() {

#ifndef USE_CUDA

  vector< vector<unsigned char> > volumeReconstructionCollection;
  volumeReconstructionCollection.push_back(mCurrentVolumeReconstruction);

  vector< vector<float> > sourceAttenuationCollection;
  sourceAttenuationCollection.push_back(mCurrentVolumeSourceAttenuation);

  vector<float> forwardProjectionAccum(1);
  vector<float> forwardProjection(1);
  vector<float> ts;
  vector< vector<unsigned char> > mats;
  vector<float> sumDetectorAttenuation;
  vector<float> startSourceAttenCollection(1);
  vector<float> endSourceAttenCollection(1);



  for (int p=0; p<mGeometry.GetNumProjectionAngles(); p++) {
    for (int y=0; y<mGeometry.GetDetectorSamplesHeight(); y++) {
      for (int x=0; x<mGeometry.GetDetectorSamplesWidth(); x++) {
        int pi = mGeometry.ProjectionCoordToIndex(p, x, y);

        mCurrentForwardProjection[pi] = 0;
        for (int s=0; s<mSamplesPerPixel; s++) {

          GetForwardProjectionCollectionForPixel(p, x, y,
                                                 (mGeometry.GetTotalProjectionSamples() * s + pi),
                                                 0, 1e5, // full ray length
                                                 1, // single volume in collection
                                                 volumeReconstructionCollection,
                                                 sourceAttenuationCollection,
                                                 forwardProjection,
                                                 ts, mats,
                                                 startSourceAttenCollection,
                                                 endSourceAttenCollection,
                                                 sumDetectorAttenuation);

        
          mCurrentForwardProjection[pi] += forwardProjection[0];
        }

        mCurrentForwardProjection[pi] /= mSamplesPerPixel;
      }
    }
  }

#else

  vector< vector<float> > forwardProjectionCollection;
  CudaForwardProject(1, &forwardProjectionCollection);
  mCurrentForwardProjection = forwardProjectionCollection[0];

#endif

}


void MarkovContext::GetAffectedProjectionPixelsParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                                        const Cone &sourceChangeCone,
                                                        vector<int> &piList) const {
  vector<int> myPiList;

  for (int pi=id; pi<mGeometry.GetTotalProjectionSamples(); pi+=nt) {
    int p, x, y;
    mGeometry.ProjectionIndexToCoord(pi, p, x, y);

    // generate all rays
    for (int s=0; s<mSamplesPerPixel; s++) {
      const Vec3f &detectorRayOrigin = mDetectorRayOrigin[mGeometry.GetTotalProjectionSamples() * s + pi];
      const Vec3f &detectorRayDir = mDetectorRayDir[mGeometry.GetTotalProjectionSamples() * s + pi];

      float t0, t1;
      if (sourceChangeCone.Intersect(detectorRayOrigin, detectorRayDir,
                                     t0, t1)) {
        myPiList.push_back(pi);
        break;
      }
    }
  }


  ParallelEnterCriticalSection(mutex);
  for (int i=0; i<myPiList.size(); i++)
    piList.push_back(myPiList[i]);
  ParallelLeaveCriticalSection(mutex);
}


void MarkovContext::UpdateForwardProjectionCollectionParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                                              volatile int workToken,
                                                              const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                                              const vector<float> &currentForwardProjection,
                                                              const Cone &sourceChangeCone,
                                                              const vector< vector<float> > &sourceAttenuationCollection,
                                                              const vector<int> &piList,
                                                              vector< vector<float> > &forwardProjectionCollection,
                                                              std::pair<int,int> &rayCasts) const {

  const int collectionSize = volumeReconstructionCollection.size();
  vector<float> forwardProjectionAccum(collectionSize);
  vector<float> forwardProjection(collectionSize);

  vector<float> t0(mSamplesPerPixel);
  vector<float> t1(mSamplesPerPixel);

  vector<float> ts;
  vector< vector<unsigned char> > mats;
  vector<float> startSourceAttenCollection(collectionSize);
  vector<float> endSourceAttenCollection(collectionSize);


  vector<float> sumDetectorAttenuation(collectionSize);

  vector<float> tforwardProjection(collectionSize);
  vector<float> tsumDetectorAttenuation(collectionSize);


  int myRayCasts = 0;
  
  for (int _pi=id; _pi<piList.size(); _pi+=nt) {
    int pi = piList[_pi];

    int p, x, y;
    mGeometry.ProjectionIndexToCoord(pi, p, x, y);

    // check if any intersect the change cone
    bool hasIntersection = false;
    for (int s=0; s<mSamplesPerPixel; s++) {
      if (sourceChangeCone.Intersect(mDetectorRayOrigin[mGeometry.GetTotalProjectionSamples() * s + pi],
                                     mDetectorRayDir[mGeometry.GetTotalProjectionSamples() * s + pi],
                                     t0[s], t1[s])) {
        hasIntersection = true;
        break;
      }
    } 

    if (hasIntersection) {
      myRayCasts++;

      for (int c=0; c<collectionSize; c++)
        forwardProjectionAccum[c] = 0;

      for (int s=0; s<mSamplesPerPixel; s++) {
#if 0
        // single cast through all collection variants
        GetForwardProjectionCollectionForPixel(p, x, y,
                                               (mGeometry.GetTotalProjectionSamples() * s + pi),
                                               0, 1e5,
                                               volumeReconstructionCollection.size(),
                                               volumeReconstructionCollection,
                                               sourceAttenuationCollection,
                                               forwardProjection,
                                               ts, mats, 
                                               startSourceAttenCollection,
                                               endSourceAttenCollection,
                                               tsumDetectorAttenuation);
        

#else
        // first cast from origin to first cone intersection, only through a single volume
        GetForwardProjectionCollectionForPixel(p, x, y,
                                               (mGeometry.GetTotalProjectionSamples() * s + pi),
                                               0, t0[s],
                                               1,
                                               volumeReconstructionCollection,
                                               sourceAttenuationCollection,
                                               tforwardProjection,
                                               ts, mats, 
                                               startSourceAttenCollection,
                                               endSourceAttenCollection,
                                               tsumDetectorAttenuation);

        for (int c=0; c<collectionSize; c++) {
          forwardProjection[c] = tforwardProjection[0];
          sumDetectorAttenuation[c] = tsumDetectorAttenuation[0];
        }
        

        // next cast through the collection variants, only inside the change area
        GetForwardProjectionCollectionForPixel(p, x, y,
                                               (mGeometry.GetTotalProjectionSamples() * s + pi),
                                               t0[s], t1[s],
                                               volumeReconstructionCollection.size(),
                                               volumeReconstructionCollection,
                                               sourceAttenuationCollection,
                                               tforwardProjection,
                                               ts, mats,
                                               startSourceAttenCollection,
                                               endSourceAttenCollection,
                                               tsumDetectorAttenuation);
        
        for (int c=0; c<collectionSize; c++) {
          forwardProjection[c] += tforwardProjection[c] * exp(sumDetectorAttenuation[c]);
          sumDetectorAttenuation[c] += tsumDetectorAttenuation[c];
        }


        // finally cast through the far side of the change, only through a single volume
        GetForwardProjectionCollectionForPixel(p, x, y,
                                               (mGeometry.GetTotalProjectionSamples() * s + pi),
                                               t1[s], 1e5,
                                               1,
                                               volumeReconstructionCollection,
                                               sourceAttenuationCollection,
                                               tforwardProjection,
                                               ts, mats, 
                                               startSourceAttenCollection,
                                               endSourceAttenCollection,
                                               tsumDetectorAttenuation);

        for (int c=0; c<collectionSize; c++) {
          forwardProjection[c] += tforwardProjection[0] * exp(sumDetectorAttenuation[c]);
        }

#endif


        for (int c=0; c<collectionSize; c++) {
          forwardProjectionAccum[c] += forwardProjection[c];
        }
      }

      for (int c=0; c<collectionSize; c++)
        forwardProjectionCollection[c][pi] = forwardProjectionAccum[c] / mSamplesPerPixel;
    }

    else {
      for (int c=0; c<collectionSize; c++)
        forwardProjectionCollection[c][pi] = currentForwardProjection[pi];
    }

  }

  ParallelEnterCriticalSection(mutex);
  rayCasts.first = std::max(rayCasts.first, myRayCasts);
  rayCasts.second += myRayCasts;
  ParallelLeaveCriticalSection(mutex);

}



// compute the error between the current forward projection and the baseline
double MarkovContext::ComputeTotalError(const vector<unsigned char> &volumeReconstruction,
                                        const vector<float> &forwardProjection) const {

  return ComputeProjectionError(forwardProjection) + ComputeRegularizationError(volumeReconstruction);
}


double MarkovContext::ComputeProjectionError(const vector<float> &forwardProjection) const {
  double projError = 0;
  for (int pi=0; pi<mGeometry.GetTotalProjectionSamples(); pi++) {
    double diff = forwardProjection[pi] - mBaselineProjection[pi];
    projError += diff*diff;
  }
  return projError * mGeometry.GetDetectorPixelArea();
}


double MarkovContext::ComputeRegularizationError(const vector<unsigned char> &volumeReconstruction) const {

  int disagreements_x = 0;
  for (int z=0; z<mGeometry.GetVolumeNodeSamplesZ(); z++) {
    for (int y=0; y<mGeometry.GetVolumeNodeSamplesY(); y++) {
      for (int x=0; x<mGeometry.GetVolumeNodeSamplesX()-1; x++) {
        if (volumeReconstruction[mGeometry.VolumeNodeCoordToIndex(x,y,z)] != 
            volumeReconstruction[mGeometry.VolumeNodeCoordToIndex(x+1,y,z)])
          disagreements_x++;
      }
    }
  }

  int disagreements_y = 0;
  for (int z=0; z<mGeometry.GetVolumeNodeSamplesZ(); z++) {
    for (int y=0; y<mGeometry.GetVolumeNodeSamplesY()-1; y++) {
      for (int x=0; x<mGeometry.GetVolumeNodeSamplesX(); x++) {
        if (volumeReconstruction[mGeometry.VolumeNodeCoordToIndex(x,y,z)] != 
            volumeReconstruction[mGeometry.VolumeNodeCoordToIndex(x,y+1,z)])
          disagreements_y++;
      }
    }
  }

  int disagreements_z = 0;
  for (int z=0; z<mGeometry.GetVolumeNodeSamplesZ()-1; z++) {
    for (int y=0; y<mGeometry.GetVolumeNodeSamplesY(); y++) {
      for (int x=0; x<mGeometry.GetVolumeNodeSamplesX(); x++) {
        if (volumeReconstruction[mGeometry.VolumeNodeCoordToIndex(x,y,z)] != 
            volumeReconstruction[mGeometry.VolumeNodeCoordToIndex(x,y,z+1)])
          disagreements_z++;
      }
    }
  }

  Vec3f voxelDim;
  mGeometry.GetVoxelDimensionsCell(voxelDim);
  return (mEnergyRegularizationWeight * (disagreements_x * voxelDim[1]*voxelDim[2] +
                                         disagreements_y * voxelDim[0]*voxelDim[2] + 
                                         disagreements_z * voxelDim[0]*voxelDim[1]));
}



// get random permutation of voxel proposals
void MarkovContext::GetGibbsIterationProposals(vector<GibbsProposal> &proposals) {

  int skipped = 0;

  proposals.reserve(mGeometry.GetTotalVolumeSamples()*3);
  proposals.clear();

  for (int z=0; z<mGeometry.GetVolumeNodeSamplesZ(); z++) {
    for (int y=0; y<mGeometry.GetVolumeNodeSamplesY(); y++) {
      for (int x=0; x<mGeometry.GetVolumeNodeSamplesX(); x++) {
        int v1 = mGeometry.VolumeNodeCoordToIndex(x,y,z);

        // use all proposals, even if they're outside right now
        /*
        bool hasSourceIntensity = false;
        for (int z2=-1; z2<2; z2+=2) {
          for (int y2=-1; y2<2; y2+=2) {
            for (int x2=-1; x2<2; x2+=2) {
              Vec3f voxelCornerPos;
              mGeometry.VolumeToWorld(Vec3f(x+x2, y+y2, z+z2), voxelCornerPos);
              if (mGeometry.GetSourceIntensityThroughPoint(0, voxelCornerPos) > 0) {
                hasSourceIntensity = true;
              }
            }
          }
        }

        if (!hasSourceIntensity) {
          skipped++;
          continue;
        }
        */


        if (x!=mGeometry.GetVolumeNodeSamplesX()-1)
          proposals.push_back(GibbsProposal(v1, mGeometry.VolumeNodeCoordToIndex(x+1,y,z)));

        if (y!=mGeometry.GetVolumeNodeSamplesY()-1)
          proposals.push_back(GibbsProposal(v1, mGeometry.VolumeNodeCoordToIndex(x,y+1,z)));

        if (z!=mGeometry.GetVolumeNodeSamplesZ()-1)
          proposals.push_back(GibbsProposal(v1, mGeometry.VolumeNodeCoordToIndex(x,y,z+1)));

      }
    }
  }

  // randomly permute the list
  for (int i=proposals.size()-1; i>=1; i--) {
    std::swap(proposals[i], proposals[rand()%(i+1)]);
  }

  std::cerr<<"keeping "<<proposals.size() << " of "<< (proposals.size()+skipped) << "total proposals"<<std::endl;
}


void MarkovContext::ProposalToVolumeCollection(const vector<unsigned char> &currentVolumeReconstruction,
                                               const GibbsProposal &proposal,
                                               vector< vector<unsigned char> > &volumeReconstructionCollection,
                                               Cone &changeCone) const {

  // single voxel change proposal
  if (proposal.second < 0) {
    int vi = proposal.first;
    int vix, viy, viz;
    mGeometry.VolumeIndexToNodeCoord(vi, vix, viy, viz);

    // get the bounding box of the change proposal
    Vec3f viMin, viMax;
    mGeometry.VolumeToWorld(Vec3f(vix-1,  viy-1,  viz-1),  viMin);
    mGeometry.VolumeToWorld(Vec3f(vix+1,  viy+1,  viz+1),  viMax);

    changeCone = Cone(mGeometry.GetSourcePosition(), viMin, viMax);


    volumeReconstructionCollection.resize(mMaterials.size());
    for (int c=0; c<mMaterials.size(); c++) {
      volumeReconstructionCollection[c] = currentVolumeReconstruction;
      volumeReconstructionCollection[c][vi ] = c;
    }
  }

  // two voxel change proposal
  else {
    int vi = proposal.first;
    int vi2 = proposal.second;
    int vix, viy, viz;
    int vi2x, vi2y, vi2z;
    mGeometry.VolumeIndexToNodeCoord(vi, vix, viy, viz);
    mGeometry.VolumeIndexToNodeCoord(vi2, vi2x, vi2y, vi2z);

    // get the bounding box of the change proposal
    Vec3f viMin, viMax;
    Vec3f vi2Min, vi2Max;
    mGeometry.VolumeToWorld(Vec3f(vix-1,  viy-1,  viz-1),  viMin);
    mGeometry.VolumeToWorld(Vec3f(vix+1,  viy+1,  viz+1),  viMax);
    mGeometry.VolumeToWorld(Vec3f(vi2x-1, vi2y-1, vi2z-1), vi2Min);
    mGeometry.VolumeToWorld(Vec3f(vi2x+1, vi2y+1, vi2z+1), vi2Max);
    Vec3f changeMin(std::min(viMin[0], vi2Min[0]),
                    std::min(viMin[1], vi2Min[1]),
                    std::min(viMin[2], vi2Min[2]));
    Vec3f changeMax(std::max(viMax[0], vi2Max[0]),
                    std::max(viMax[1], vi2Max[1]),
                    std::max(viMax[2], vi2Max[2]));

    changeCone = Cone(mGeometry.GetSourcePosition(), changeMin, changeMax);


    volumeReconstructionCollection.resize(mMaterials.size()*mMaterials.size());
    for (int c=0; c<mMaterials.size()*mMaterials.size(); c++) {
      volumeReconstructionCollection[c] = currentVolumeReconstruction;
      volumeReconstructionCollection[c][vi ] = c % mMaterials.size();
      volumeReconstructionCollection[c][vi2] = c / mMaterials.size();
    }
  }
}


int MarkovContext::GibbsEval(ParallelThreadPool *threadPool, int numThreads, float temp, const GibbsProposal &proposal,
                             const vector<unsigned char> &currentVolumeReconstruction,
                             const vector<float> &currentVolumeSourceAttenuation,
                             const vector<float> &currentForwardProjection,
                             vector< vector<float> > &sourceAttenuationCollection,
                             vector< vector<float> > &forwardProjectionCollection,
                             vector<double> &matErrors,
                             vector<double> &matProbs) const {

  // get all of the different material configurations we will be considering
  Cone changeCone;
  vector< vector<unsigned char> > volumeReconstructionCollection;
  ProposalToVolumeCollection(currentVolumeReconstruction, proposal,
                             volumeReconstructionCollection, changeCone);

  int numConfigs = volumeReconstructionCollection.size();

#ifndef USE_CUDA
  // compute source attenuation of each potential change
  Cone attenChangeCone = changeCone;
  if (threadPool) {
    threadPool->ParallelRun(makeClassFunctor(this, &MarkovContext::UpdateSourceAttenuationCollectionParallel),
                            volumeReconstructionCollection, currentVolumeSourceAttenuation,
                            changeCone, attenChangeCone, sourceAttenuationCollection);
  }
  else {
    ParallelExecutor(numThreads,
                     makeClassFunctor(this, &MarkovContext::UpdateSourceAttenuationCollectionParallel),
                     volumeReconstructionCollection, currentVolumeSourceAttenuation,
                     changeCone, attenChangeCone, sourceAttenuationCollection);
  }


  // compute forward projections of each potential change
  forwardProjectionCollection.resize(numConfigs);
  for (int c=0; c<numConfigs; c++) {
    forwardProjectionCollection[c] = currentForwardProjection;
  }


  std::pair<int,int> rayCasts(0,0);

  if (threadPool) {
    vector<int> piList;
    threadPool->ParallelRun(makeClassFunctor(this, &MarkovContext::GetAffectedProjectionPixelsParallel),
                            attenChangeCone, piList);


    volatile int workToken = 0;
    threadPool->ParallelRun(makeClassFunctor(this, &MarkovContext::UpdateForwardProjectionCollectionParallel),
                            workToken,
                            volumeReconstructionCollection,
                            currentForwardProjection,
                            attenChangeCone, 
                            sourceAttenuationCollection,
                            piList,
                            forwardProjectionCollection,
                            rayCasts);
  }
  else {
    vector<int> piList;
    ParallelExecutor(numThreads,
                     makeClassFunctor(this, &MarkovContext::GetAffectedProjectionPixelsParallel),
                     attenChangeCone, piList);

    volatile int workToken = 0;
    ParallelExecutor(numThreads, 
                     makeClassFunctor(this, &MarkovContext::UpdateForwardProjectionCollectionParallel),
                     workToken,
                     volumeReconstructionCollection,
                     currentForwardProjection,
                     attenChangeCone,
                     sourceAttenuationCollection,
                     piList,
                     forwardProjectionCollection,
                     rayCasts);
  }


  // compute errors of each potential change
  matErrors.resize(numConfigs);
  for (int c=0; c<numConfigs; c++) {
    matErrors[c] = ComputeTotalError(volumeReconstructionCollection[c],
                                     forwardProjectionCollection[c]);
  }


#else

  // ray cast in cuda
  //CudaSetCurrentVolume(currentVolumeReconstruction);
  CudaSetVolumeCollection(proposal);

  /*
  // test that we set the volumes correctly
  vector< vector<unsigned char> > volCollection;
  CudaGetVolumeCollection(volCollection);
  for (int c=0; c<volCollection.size(); c++) {
    for (int i=0; i<volCollection[c].size(); i++) {
      if (volCollection[c][i] != volumeReconstructionCollection[c][i]) {
        std::cerr<<"volume collection not set correctly!"<<std::endl;
        std::cerr<<c<<" "<<i<<std::endl;
      }
    }
  }
  */


  CudaComputeSourceAttenuation(numConfigs, NULL);



  Cone attenChangeCone = changeCone;
  GetAttenChangeCone(changeCone, attenChangeCone);
  

  // compute errors and do reduction on gpu (prevents copying forward projections back to host memory)
  /*
  //CudaForwardProject(numConfigs, NULL);
  CudaUpdateForwardProjection(numConfigs, attenChangeCone, NULL);
  vector<float> cudaErrors;
  CudaGetProjectionError(numConfigs, cudaErrors);
  matErrors.resize(numConfigs);
  for (int c=0; c<numConfigs; c++) {
    matErrors[c] = cudaErrors[c];// + ComputeRegularizationError(volumeReconstructionCollection[c]);
  }
  */

  // copy projections back to host and compute errors on cpu
  CudaForwardProject(numConfigs, &forwardProjectionCollection);
  //CudaUpdateForwardProjection(numConfigs, attenChangeCone, &forwardProjectionCollection);

  matErrors.resize(numConfigs);
  for (int c=0; c<numConfigs; c++) {
    //    matErrors[c] = ComputeTotalError(volumeReconstructionCollection[c],
    //                                     forwardProjectionCollection[c]);
    matErrors[c] = ComputeProjectionError(forwardProjectionCollection[c])
      + ComputeRegularizationError(volumeReconstructionCollection[c]);
  }

  /*
  static int counter=0;
  counter++;
  if (counter==1000) {
    CudaShutdown();
    exit(0);
  }
  */
  //return 0;



  // test update vs complete recompute
  /*
  CudaForwardProject(numConfigs, &forwardProjectionCollection);
  vector< vector<float> > uforwardProjectionCollection;
  CudaUpdateForwardProjection(numConfigs, currentForwardProjection, attenChangeCone, &uforwardProjectionCollection);
  for (int c=0; c<numConfigs; c++) {
    for (int i=0; i<mGeometry.GetTotalProjectionSamples(); i++) {
      if (forwardProjectionCollection[c][i] != uforwardProjectionCollection[c][i]) {
        std::cerr<<forwardProjectionCollection[c][i]<<" "<<uforwardProjectionCollection[c][i]<<std::endl;
      }
    }
  }
  */


#endif


  // normalize probabilities, convert to cumulative probabilites
  bool print = false;
  if (print) {
    std::cerr<<"orig: ";
    for (int i=0; i<(int)matErrors.size(); i++) {
      std::cerr<<matErrors[i]<<" ";
    }
    std::cerr<<std::endl;
  }

  double minError = 1e99;
  double maxError = 0;
  for (int i=0; i<(int)matErrors.size(); i++) {
    minError = std::min(minError, matErrors[i]);
    maxError = std::max(maxError, matErrors[i]);
  }


  vector<double> matErrorsAdj(numConfigs);
  for (int i=0; i<(int)matErrors.size(); i++) {
    matErrorsAdj[i] = (matErrors[i]-minError);
  }

  //std::cerr<<"min/max error ratio: "<<minError<<" "<<maxError<<" "<<(minError/maxError)<<std::endl;
  if (print) {

    std::cerr<<"maxed: ";
    for (int i=0; i<(int)matErrors.size(); i++) {
      std::cerr<<matErrorsAdj[i]<<" ";
    }
    std::cerr<<std::endl;

    std::cerr<<"perror: ";
    for (int i=0; i<(int)matErrors.size(); i++) {
      std::cerr<<(matErrorsAdj[i]/temp)<<" ";
    }
    std::cerr<<std::endl;
  }


  // convert errors to probabilities
  matProbs.resize(numConfigs);
  for (int i=0; i<(int)matErrors.size(); i++) {
    matProbs[i] = exp(-(matErrorsAdj[i])/temp);
  }


  if (print) {
    std::cerr<<"probs: ";
    for (int i=0; i<(int)matErrors.size(); i++) {
      std::cerr<<matProbs[i]<<" ";
    }
    std::cerr<<std::endl;
  }


  double totalProb = 0;
  for (int i=0; i<(int)matErrors.size(); i++) {
    totalProb += matProbs[i];
  }
  for (int i=0; i<(int)matErrors.size(); i++) {
    matProbs[i] /= totalProb;
  }

  if (print) {
    std::cerr<<"nprobs: ";
    for (int i=0; i<(int)matErrors.size(); i++) {
      std::cerr<<matProbs[i]<<" ";
    }
    std::cerr<<std::endl;
    std::cerr<<std::endl;
  }


  // sample proportional to probabilities
  double r = (double)rand() / RAND_MAX;

  int nextConfig = 0;
  double probSum=0;
  for (int i=0; i<(int)matErrors.size(); i++) {
    probSum += matProbs[i];
    if (r < probSum) {
      nextConfig = i;
      break;
    }
  }

  return nextConfig;
}


void MarkovContext::GibbsIterParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                      volatile int &syncToken, volatile int &workToken,
                                      ParallelThreadPool *threadPool, int numSubThreads, float temp,
                                      const vector<GibbsProposal> &proposals) {

  int lastUpdate;
  ParallelEnterCriticalSection(mutex);
  lastUpdate = syncToken;
  vector<unsigned char> currentVolumeReconstruction = mCurrentVolumeReconstruction;
  vector<float> currentVolumeSourceAttenuation = mCurrentVolumeSourceAttenuation;
  vector<float> currentForwardProjection = mCurrentForwardProjection;
  ParallelLeaveCriticalSection(mutex);


  vector< vector<float> > sourceAttenuationCollection(mMaterials.size()*mMaterials.size());
  vector< vector<float> > forwardProjectionCollection(mMaterials.size()*mMaterials.size());
  vector<double> matErrors(mMaterials.size()*mMaterials.size());
  vector<double> matProbs(mMaterials.size()*mMaterials.size());

  for (int i=0; i<mMaterials.size()*mMaterials.size(); i++) {
    sourceAttenuationCollection[i].resize(mGeometry.GetTotalVolumeNodeSamples());
    forwardProjectionCollection[i].resize(mGeometry.GetTotalProjectionSamples());
  }

  while (1) {
    int _vi = -1;
    ParallelEnterCriticalSection(mutex);
    if (workToken < proposals.size()) {
      _vi = workToken;
      workToken++;
    }
    ParallelLeaveCriticalSection(mutex);

    if (_vi < 0)
      break;

    int vi = proposals[_vi].first;
    int vi2 = proposals[_vi].second;

    int vix, viy, viz;
    mGeometry.VolumeIndexToNodeCoord(vi, vix, viy, viz);

    int vi2x, vi2y, vi2z;
    mGeometry.VolumeIndexToNodeCoord(vi2, vi2x, vi2y, vi2z);

    while (1) {
      if (syncToken != lastUpdate) {
        ParallelEnterCriticalSection(mutex);
        lastUpdate = syncToken;
        currentVolumeReconstruction = mCurrentVolumeReconstruction;
        currentVolumeSourceAttenuation = mCurrentVolumeSourceAttenuation;
        currentForwardProjection = mCurrentForwardProjection;
        ParallelLeaveCriticalSection(mutex);
      }


      int currentConfig = currentVolumeReconstruction[vi] + mMaterials.size()*currentVolumeReconstruction[vi2];

      int nextConfig = GibbsEval(threadPool ? &threadPool[id] : NULL, numSubThreads, temp, proposals[_vi],
                                 currentVolumeReconstruction,
                                 currentVolumeSourceAttenuation,
                                 currentForwardProjection,
                                 sourceAttenuationCollection, forwardProjectionCollection,
                                 matErrors, matProbs);


      std::cerr<<id<<" "<<_vi<<" "<<vi<<" "<<vi2<<" "<<temp<<" "<<matErrors[nextConfig]<<
        " "<<currentConfig<<" "<<nextConfig<<std::endl;

      currentVolumeReconstruction[vi] = nextConfig % mMaterials.size();
      currentVolumeReconstruction[vi2] = nextConfig / mMaterials.size();
      currentVolumeSourceAttenuation = sourceAttenuationCollection[nextConfig];
      currentForwardProjection = forwardProjectionCollection[nextConfig];

#ifdef USE_CUDA
      CudaAcceptNextConfig(nextConfig);
#endif


      bool proposalDone = true;
      if (nextConfig != currentConfig) {
        ParallelEnterCriticalSection(mutex);
        if (lastUpdate == syncToken) {
          syncToken++;
          lastUpdate = syncToken;
          mCurrentVolumeReconstruction = currentVolumeReconstruction;
          mCurrentVolumeSourceAttenuation = currentVolumeSourceAttenuation;
          mCurrentForwardProjection = currentForwardProjection;
        }
        else {
          // if we failed to update, try the same change with new base configuration
          proposalDone = false;
        }
        ParallelLeaveCriticalSection(mutex);
      }

      if (proposalDone)
        break;
    }
  }
}


// run Gibbs sampler
void MarkovContext::Gibbs(int numThreads, int numSubThreads, int niter, float startTemp, float endTemp) {

#ifdef USE_CUDA
  numThreads = 1;
  numSubThreads = 1;
#endif

  // get the current energy
  ComputeSourceAttenuation();
  ComputeForwardProjection();

  ParallelThreadPool *threadPool = new ParallelThreadPool[numThreads];
  for (int i=0; i<numThreads; i++) {
    threadPool[i].StartThreads(numSubThreads);
  }
  
  for (int iter=0; iter<niter; iter++) {

    double T = startTemp * (niter-iter)/niter + endTemp * iter/niter;
    T = startTemp / log(endTemp*iter+2.0);

    // create a voxel permutation
    vector<GibbsProposal> proposals;
    GetGibbsIterationProposals(proposals);

    volatile int syncToken=0, workToken=0;
    ParallelExecutor(numThreads, 
                     makeClassFunctor(this, &MarkovContext::GibbsIterParallel),
                     syncToken, workToken, threadPool, numSubThreads, T, proposals);
  }

  for (int i=0; i<numThreads; i++) {
    threadPool[i].StopThreads();
  }
}



void MarkovContext::GetGibbsConditionalProbabilities(int numThreads, float temp, 
                                                     vector< vector<float> > &errors,
                                                     vector< vector<float> > &probs) const {
  ParallelThreadPool threadPool;
  threadPool.StartThreads(numThreads);


  std::cerr<<"computing reconstruction probabilities"<<std::endl;

  errors.resize(mMaterials.size());
  probs.resize(mMaterials.size());
  for (int i=0; i<probs.size(); i++) {
    errors[i].resize(mCurrentVolumeReconstruction.size());
    probs[i].resize(mCurrentVolumeReconstruction.size());
  }


  vector< vector<float> > sourceAttenuationCollection(mMaterials.size()*mMaterials.size());
  vector< vector<float> > forwardProjectionCollection(mMaterials.size()*mMaterials.size());
  vector<double> matErrors(mMaterials.size()*mMaterials.size());
  vector<double> matProbs(mMaterials.size()*mMaterials.size());

  for (int i=0; i<mMaterials.size()*mMaterials.size(); i++) {
    sourceAttenuationCollection[i].resize(mGeometry.GetTotalVolumeNodeSamples());
    forwardProjectionCollection[i].resize(mGeometry.GetTotalProjectionSamples());
  }


  for (int vi=0; vi<mCurrentVolumeReconstruction.size(); vi++) {

    //std::cerr<<"reconstruction probabilities: "<<vi<<std::endl;

    /*
    GibbsProposal proposal(vi, -1);
    GibbsEval(&threadPool, 0, temp, proposal,
              mCurrentVolumeReconstruction,
              mCurrentVolumeSourceAttenuation,
              mCurrentForwardProjection,
              sourceAttenuationCollection,
              forwardProjectionCollection,
              matErrors,
              matProbs);
    */
    matErrors = std::vector<double>(matErrors.size(), 1);
    matProbs = std::vector<double>(matProbs.size(), 1);



    float tprob = 0;
    for (int m=0; m<mMaterials.size(); m++)
      tprob += matProbs[m];
    for (int m=0; m<mMaterials.size(); m++)
      probs[m][vi] = matProbs[m] / tprob;
    for (int m=0; m<mMaterials.size(); m++)
      errors[m][vi] = matErrors[m];
  }



  threadPool.StopThreads();
}

