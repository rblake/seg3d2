
#ifndef MARKOV_H
#define MARKOV_H

#define NUM_MATERIALS 3

#include <vector>

#include <Application/BackscatterReconstruction/Algorithm/vec3.h>
#include <Application/BackscatterReconstruction/Algorithm/cone.h>
#include <Application/BackscatterReconstruction/Algorithm/geometry.h>
#include <Application/BackscatterReconstruction/Algorithm/material.h>


#ifndef __CUDACC__
#include <Application/BackscatterReconstruction/Algorithm/parallel.h>
#include <itkImage.h>
typedef itk::Image<float, 3> FloatVolumeType;
typedef itk::Image<unsigned char, 3> ByteVolumeType;
#endif


typedef std::pair<int,int> GibbsProposal;

class MarkovContext {
public:

#ifndef __CUDACC__
  MarkovContext(const Geometry &g, const vector<Material> &materials, 
                int samplesPerPixel, float voxelStepSize, float energyRegularizationWeight);

  // set current volume reconstruction
  void SetCurrentVolume(ByteVolumeType::Pointer vol);
  void SetCurrentVolume(unsigned char v);

  // get current volume reconstruction
  void GetCurrentVolume(ByteVolumeType::Pointer vol);

  // set baseline projections used for reconstruction/backward projection
  void SetBaselineProjection(FloatVolumeType::Pointer vol);

  // get current projection
  void GetCurrentProjection(FloatVolumeType::Pointer vol);

  // get current source attenuation
  void GetCurrentVolumeSourceAttenuation(FloatVolumeType::Pointer vol) const;

  // run simple optimization to find a decent initial guess
  void FindInitialGuess();

  // integrate the attenuation coefficients from the source into the volume
  void ComputeSourceAttenuation();

  void UpdateSourceAttenuationCollectionParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                                 const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                                 const vector<float> &currentVolumeSourceAttenuation,
                                                 const Cone &sourceChangeCone,
                                                 Cone &attenChangeCone,
                                                 vector< vector<float> > &sourceAttenuation) const;

  void UpdateSourceAttenuationCollectionSample(const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                               const vector<float> &currentVolumeSourceAttenuation,
                                               const Cone *sourceChangeCone,
                                               int nvi,
                                               Cone &attenChangeCone,
                                               vector< vector<float> > &sourceAttenuation) const;


  // use the current volume guess and source attenuation to do a forward projection
  void ComputeForwardProjection();

  void GetAffectedProjectionPixelsParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                           const Cone &sourceChangeCone,
                                           vector<int> &piList) const;

  void UpdateForwardProjectionCollectionParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                                                 volatile int workToken,
                                                 const vector< vector<unsigned char> > &volumeReconstructionCollection,
                                                 const vector<float> &currentForwardProjection,
                                                 const Cone &sourceChangeCone,
                                                 const vector< vector<float> > &sourceAttenuationCollection,
                                                 const vector<int> &piList,
                                                 vector< vector<float> > &forwardProjectionCollection,
                                                 std::pair<int,int> &rayCasts) const;

  void GetForwardProjectionCollectionForPixel(int p, int x, int y,
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
                                              vector<float> &sumDetectorAttenuation) const;


  // run Gibbs sampler
  int GibbsEval(ParallelThreadPool *threadPool, int numThreads, float temp, 
                const GibbsProposal &proposal, const Vec3i &baseRegularizationCounts,
                const vector<unsigned char> &currentVolumeReconstruction,
                const vector<float> &currentVolumeSourceAttenuation,
                const vector<float> &currentForwardProjection,
                vector< vector<float> > &sourceAttenuationCollection,
                vector< vector<float> > &forwardProjectionCollection,
                vector<double> &matErrors,
                vector<double> &matProbs) const;
  void Gibbs(int numThreads, int numSubThreads, int niter, float startTemp, float endTemp);
  void GibbsIterParallel(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex,
                         volatile int &syncToken, volatile int &workToken,
                         ParallelThreadPool *threadPool, int numSubThreads, float temp,
                         const vector<GibbsProposal> &proposals);
  void GetGibbsIterationProposals(vector<GibbsProposal> &proposals);

  void GetGibbsConditionalProbabilities(int numThreads, float temp, 
                                        vector< vector<float> > &errors,
                                        vector< vector<float> > &probs) const;
#endif

private:

#ifndef __CUDACC__
  // compute the error between the current forward projection and the baseline
  double ComputeTotalError(const vector<unsigned char> &volumeReconstruction,
                           const vector<float> &forwardProjection) const;
  double ComputeProjectionError(const vector<float> &forwardProjection) const;

  void ComputeFullRegularizationCounts(const vector<unsigned char> &volumeReconstruction, Vec3i &counts) const;
  double ComputeFullRegularizationError(const vector<unsigned char> &volumeReconstruction) const;
  double ComputeRegularizationErrorFromCounts(const Vec3i &counts) const;
  void UpdateRegularizationNeighborCounts(const vector<unsigned char> &volumeReconstruction,
                                          const GibbsProposal &proposal, int c,
                                          const Vec3i &prevCounts,
                                          Vec3i &newCounts) const;
                                          

  // get a list of intervals corresponding to marching the given ray through the volume
  bool GetRayVoxelIntersections(const vector< vector<unsigned char> > &volumeReconstructionCollection, 
                                int collectionSize,
                                int rayIndex, bool detectorRay,
                                float tmin, float tmax,
                                vector<float> &ts, vector< vector<unsigned char> > &materials) const;

  void LerpNodeValuesCollection(const vector< vector<float> > &valueCollection,
                                int collectionSize,
                                const Vec3f &volumePos,
                                vector<float> &lerpedValues) const;

  void ProposalToVolumeCollection(const vector<unsigned char> &currentVolumeReconstruction,
                                  const GibbsProposal &proposal,
                                  vector< vector<unsigned char> > &volumeReconstructionCollection,
                                  Cone &changeCone) const;

  void GetAttenChangeCone(const Cone &changeCone, Cone &attenChangeCone) const;
  
#endif

  void CudaInitialize();
  void CudaShutdown() const;
  void CudaSetBaselineProjection() const;

  void CudaSetCurrentVolume(const vector<unsigned char> &matids) const;
  void CudaSetVolumeCollection(const GibbsProposal &proposal) const;
  void CudaGetVolumeCollection(vector< vector<unsigned char> > &volumeCollection) const;
  void CudaAcceptNextConfig(int c) const;

  void CudaComputeSourceAttenuation(int collectionSize) const;
  void CudaGetSourceAttenuation(int collectionSize,
                                vector< vector<float> > &sourceAttenuationCollection) const;

  void CudaForwardProject(int collectionSize) const;
  void CudaUpdateForwardProjection(int collectionSize,
                                   const Cone &attenChangeCone) const;
  void CudaGetForwardProjection(int collectionSize,
                                vector< vector<float> > &forwardProjectionCollection) const;

  void CudaGetProjectionError(int collectionSize, vector<float> &errors) const;


  // the geometric configuration of the device
  const Geometry &mGeometry;

  const int mSamplesPerPixel;
  const float mVoxelStepSize;

  float mProgressIncrement;

  // energy regularization weight
  const float mEnergyRegularizationWeight;

  // the material attributes for all materials in the volume
  const vector<Material> &mMaterials;

  // the captured projections that we want to match
  vector<float> mBaselineProjection;

  // the current forward projection
  vector<float> mCurrentForwardProjection;

  // the current guess for the volume state
  vector<unsigned char> mCurrentVolumeReconstruction;

  // how much source energy reaches each voxel in the volume given the current volume state
  vector<float> mCurrentVolumeSourceAttenuation;

  // generate fixed set of rays to use
  void InitRay(const Vec3f &origin, const Vec3f &dir,
               Vec3f &volOrigin, Vec3f &volDir,
               float &tmin, float &tmax) const;
  vector<Vec3f> mDetectorRayOrigin;
  vector<Vec3f> mDetectorRayDir;
  vector<Vec3f> mDetectorRayVolOrigin;
  vector<Vec3f> mDetectorRayVolDir;
  vector<float> mDetectorRayTMin;
  vector<float> mDetectorRayTMax;

  vector<Vec3f> mSourceRayOrigin;
  vector<Vec3f> mSourceRayDir;
  vector<Vec3f> mSourceRayVolOrigin;
  vector<Vec3f> mSourceRayVolDir;
  vector<float> mSourceRayTMin;
  vector<float> mSourceRayTMax;
};


#endif

