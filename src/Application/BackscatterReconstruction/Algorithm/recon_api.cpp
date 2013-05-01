
#include "markov.h"
#include "calibration.h"
#include "recon_api.h"


// calibration - segment image stack of calibration pattern to find disks
void CalibrationSegment(const ReconImageVolumeType::Pointer images,
                        ReconMaterialIdVolumeType::Pointer diskIds) {


  ReconImageVolumeType::SizeType inSize = images->GetLargestPossibleRegion().GetSize();

  // create output volume
  ReconMaterialIdVolumeType::SizeType size;
  size.SetElement(0, inSize[0]);
  size.SetElement(1, inSize[1]);
  size.SetElement(2, inSize[2]);
  diskIds->SetRegions(ByteVolumeType::RegionType(size));
  diskIds->Allocate();


  // copy input into vector
  vector<float> forward(inSize[0]*inSize[1]*inSize[2]);
  for (size_t i=0; i<forward.size(); i++) {
    forward[i] = images->GetBufferPointer()[i];
  }


  // do the segmentation
  vector<unsigned char> ids;
  int numPoints = 3;
  int maskSize = 14 * inSize[0] / 64;
  int idSize = 6 * inSize[0] / 64;
  int border = 10 * inSize[0] / 64;
  bool dark = false;
  SegmentCalibPoints(inSize[0], inSize[1], inSize[2],
                     forward, ids, numPoints, maskSize, idSize, border, dark);


  // copy result back
  for (size_t i=0; i<ids.size(); i++) {
    diskIds->GetBufferPointer()[i] = ids[i];
  }
}


// calibration - given images and disk id masks (after mapping to match calibration pattern), fit detector geometry
void CalibrationFitGeometry(const ReconImageVolumeType::Pointer images,
                            const ReconMaterialIdVolumeType::Pointer diskIds,
                            const char *workDirectory,
                            const char *initialGeometryConfigFile,
                            const char *outputGeometryConfigFile) {

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(initialGeometryConfigFile);


  //
  // localize the disks
  //

  ReconImageVolumeType::SizeType inSize = images->GetLargestPossibleRegion().GetSize();

  // copy input into vectors
  vector<float> forward(inSize[0]*inSize[1]*inSize[2]);
  for (size_t i=0; i<forward.size(); i++) {
    forward[i] = images->GetBufferPointer()[i];
  }

  vector<unsigned char> maskv(inSize[0]*inSize[1]*inSize[2]);
  for (size_t i=0; i<maskv.size(); i++) {
    maskv[i] = diskIds->GetBufferPointer()[i];
  }


  // run localization
  vector<int> proj_i, point_i;
  vector<float> point_x, point_y;
  vector<float> residuals;
  LocalizeCalibPoints(inSize[0], inSize[1], inSize[2],
                      forward, maskv,
                      proj_i, point_i, point_x, point_y,
                      residuals,
                      false);

  std::cerr<<"found "<<proj_i.size()<<"points"<<std::endl;


  // 
  // fit the detector position
  //

  vector<Vec3f> cpoints;
  cpoints.push_back(Vec3f( 0.03175,-0.03175, 0));
  cpoints.push_back(Vec3f(-0.03175, 0.03175, 0));
  cpoints.push_back(Vec3f(-0.03175,-0.03175, 0));


  // do the linear algebra
  Vec3f baseCenter, baseX, baseY, baseZ;
  CalibrateDetector(geometry, proj_i, point_i, point_x, point_y, cpoints,
                    baseCenter, baseX, baseY, baseZ);


  // save the new geometry file
  geometry.SaveToFile(outputGeometryConfigFile);

}


// reconstruction - don't return until reconstruction is completed
void ReconstructionStart(const ReconImageVolumeType::Pointer images,
                         const ReconMaterialIdVolumeType::Pointer initialGuess, // possibly NULL
                         const char *geometryConfigFile,
                         int iterations) {

}



// reconstruction - stop it if it's currently running
void ReconstructionAbort() {

}


// reconstruction - get the latest reconstructed volume material ids, can be called any time during reconstruction
void ReconstructionGetMaterialVolume(ReconMaterialIdVolumeType::Pointer reconVolume) {

}

