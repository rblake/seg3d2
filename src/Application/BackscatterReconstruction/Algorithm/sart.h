#ifndef SART_H
#define SART_H

#include <vector>
#include <utility>

#include <itkImage.h>

#include <Application/BackscatterReconstruction/Algorithm/vec3.h>
#include <Application/BackscatterReconstruction/Algorithm/geometry.h>
#include <Application/BackscatterReconstruction/Algorithm/parallel.h>
#include <Application/BackscatterReconstruction/Algorithm/material.h>

using std::vector;
using std::pair;

typedef itk::Image<float, 3> FloatVolumeType;
typedef itk::Image<unsigned char, 3> ByteVolumeType;


class SARTContext {
public:
  SARTContext(const Geometry &g, const Material &foam) : mGeometry(g), mFoamMaterial(foam) {
    mBaselineProjection.resize(g.GetTotalProjectionSamples(), 0);
    mCurrentForwardProjection.resize(g.GetTotalProjectionSamples(), 0);
    mCurrentForwardProjectionDeltaR.resize(g.GetTotalProjectionSamples(), 0);
    mCurrentVolumeReconstruction.resize(g.GetTotalVolumeSamples(), 0);
  }

  // set current volume reconstruction
  void SetCurrentVolume(FloatVolumeType::Pointer vol);
  void SetCurrentVolume(float v);

  // set baseline projections used for reconstruction/backward projection
  void SetBaselineProjection(FloatVolumeType::Pointer vol);

  // get current projection
  void GetCurrentProjection(FloatVolumeType::Pointer vol);

  // get current delta r
  void GetCurrentDeltaR(FloatVolumeType::Pointer vol);

  // get current volume reconstruction
  void GetCurrentVolume(FloatVolumeType::Pointer vol);

  // compute current weights for the volume
  void CreateWeights(float stepSize, int samplesPerPixel, int numThreads);
  void CreateWeightsParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                             float stepSize, int samplesPerPixel);

  // use the weights / volume to create a forward projection
  void ForwardProject();

  // use the weights / projecitons to create a forward projection
  void BackwardProject(float lambda, int samplesPerPixel);

private:

  // the geometric configuration of the device
  const Geometry &mGeometry;

  // the material attributes for the foam
  const Material &mFoamMaterial;

  // the captured projections that we want to match
  vector<float> mBaselineProjection;

  // the current forward projection
  vector<float> mCurrentForwardProjection;
  vector<float> mCurrentForwardProjectionDeltaR;

  // the current backward projection
  vector<float> mCurrentVolumeReconstruction;


  // projection/volume weight pairs - stored mapping both ways
  vector< vector< pair<int,float> > > mVolumeToProjectionWeights;
  vector< vector< pair<int,float> > > mProjectionToVolumeWeights;


};



#endif
