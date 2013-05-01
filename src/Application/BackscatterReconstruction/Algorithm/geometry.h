#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <fstream>
#include <vector>

#include "vec3.h"

using std::vector;

class Geometry {

public:
  Geometry() { }

  bool LoadFromFile(const char *fname);
  bool SaveToFile(const char *fname) const;


  // convert projection index types
  int ProjectionCoordToIndex(int p, int x, int y) const;
  void ProjectionIndexToCoord(int i, int &p, int &x, int &y) const;
  int GetProjectionStrideX() const { return mProjectionStride[0]; }
  int GetProjectionStrideY() const { return mProjectionStride[1]; }
  int GetProjectionStrideP() const { return mProjectionStride[2]; }
  int GetProjectionStride(int axis) const { return mProjectionStride[axis]; }

  void ProjectionPixelToPhysicalOffset(float px, float py, float &mx, float &my) const;
  void ProjectionPhysicalOffsetToPixel(float mx, float my, float &px, float &py) const;


  // convert volume index types
  int VolumeCoordToIndex(int x, int y, int z) const;
  void VolumeIndexToCoord(int i, int &x, int &y, int &z) const;
  int GetVolumeStrideX() const { return mVolumeStride[0]; }
  int GetVolumeStrideY() const { return mVolumeStride[1]; }
  int GetVolumeStrideZ() const { return mVolumeStride[2]; }
  int GetVolumeStride(int axis) const { return mVolumeStride[axis]; }

  int GetVolumeNodeStrideX() const { return mVolumeNodeStride[0]; }
  int GetVolumeNodeStrideY() const { return mVolumeNodeStride[1]; }
  int GetVolumeNodeStrideZ() const { return mVolumeNodeStride[2]; }
  int GetVolumeNodeStride(int axis) const { return mVolumeNodeStride[axis]; }

  int VolumeNodeCoordToIndex(int x, int y, int z) const;
  void VolumeIndexToNodeCoord(int i, int &x, int &y, int &z) const;


  // map point from world coordinates to volume pixel coordinates, return false if outside the volume.
  // map point from volume pixel coordinates to world coordinates
  bool WorldToVolume(const Vec3f &world, Vec3f &volume) const;
  void VolumeToWorld(const Vec3f &volume, Vec3f &world) const;
  void WorldToVolumeDir(const Vec3f &worldDir, Vec3f &volumeDir) const;
  void GetVoxelDimensionsCell(Vec3f &dim) const { dim = Vec3f((mVolumeBounds[1][0]-mVolumeBounds[0][0]) / mVolumeSamples[0],
                                                              (mVolumeBounds[1][1]-mVolumeBounds[0][1]) / mVolumeSamples[1],
                                                              (mVolumeBounds[1][2]-mVolumeBounds[0][2]) / mVolumeSamples[2]); }
  void GetVoxelDimensionsNode(Vec3f &dim) const { dim = Vec3f((mVolumeBounds[1][0]-mVolumeBounds[0][0]) / (mVolumeSamples[0]-1),
                                                              (mVolumeBounds[1][1]-mVolumeBounds[0][1]) / (mVolumeSamples[1]-1),
                                                              (mVolumeBounds[1][2]-mVolumeBounds[0][2]) / (mVolumeSamples[2]-1)); }


  // return the min and max distance along the ray where it intersects with the volume's bounding box.
  // tmin/tmax may be negative!  Return false if no intersection.
  bool RayVolumeIntersection(const Vec3f &origin, const Vec3f &dir,
                             float &tmin, float &tmax) const;


  // for a point in world space, find the voxels and weights for trilinear sampling at that point.
  // return false if outside of the volume
  bool WorldToTrilinearWeights(const Vec3f &world, int indices[8], float weights[8]) const;


  int GetNumProjectionAngles() const { return (int)mProjectionAngles.size(); }
  float GetProjectionAngle(int p) const { return mProjectionAngles[p]; }
  int GetDetectorSamplesWidth() const { return mDetectorSamples[0]; }
  int GetDetectorSamplesHeight() const { return mDetectorSamples[1]; }
  int GetTotalProjectionSamples() const { return mDetectorSamples[0] * mDetectorSamples[1] * (int)mProjectionAngles.size(); }

  int GetVolumeSamplesX() const { return mVolumeSamples[0]; }
  int GetVolumeSamplesY() const { return mVolumeSamples[1]; }
  int GetVolumeSamplesZ() const { return mVolumeSamples[2]; }
  int GetVolumeNodeSamplesX() const { return mVolumeSamples[0]+1; }
  int GetVolumeNodeSamplesY() const { return mVolumeSamples[1]+1; }
  int GetVolumeNodeSamplesZ() const { return mVolumeSamples[2]+1; }
  int GetTotalVolumeSamples() const { return mVolumeSamples[0] * mVolumeSamples[1] * mVolumeSamples[2]; }
  int GetTotalVolumeNodeSamples() const { return (mVolumeSamples[0]+1) * (mVolumeSamples[1]+1) * (mVolumeSamples[2]+1); }

  Vec3f GetSourcePosition() const;
  Vec3f GetSourceDirection() const;
  float GetSourceIntensityThroughPoint(const Vec3f &pos) const;

  // generate a random ray that hits the pixel, and passes through the entire collimation column
  float GetColimatorFocalLength() const { return mDetectorCollimatorFocalLength; }
  void GetDetectorRay(int p, int x, int y, bool centered,
                      Vec3f &origin, Vec3f &dir) const;

  float GetDetectorPixelArea() const { return ((mDetectorSize[0]/ mDetectorSamples[0]) *
                                               (mDetectorSize[1]/ mDetectorSamples[1])); }


  template <typename T>
  Vec3<T> GetRotatedPoint(const Vec3<T> &vec, int i) const { return RotateVector(vec, mProjectionAngles[i]); }

  template <typename T>
  Vec3<T> GetInverseRotatedPoint(const Vec3<T> &vec, int i) const { return RotateVector(vec, -mProjectionAngles[i]); }


  Vec3f GetDetectorForward() const { return mDetectorForward; }
  Vec3f GetDetectorRight() const { return mDetectorRight; }
  Vec3f GetDetectorUp() const { return mDetectorUp; }
  Vec3f GetDetectorCenter() const { return mDetectorCenter; }

  void SetDetectorForward(const Vec3f &x) { mDetectorForward = x; }
  void SetDetectorRight(const Vec3f &x) { mDetectorRight = x; }
  void SetDetectorUp(const Vec3f &x) { mDetectorUp = x; }
  void SetDetectorCenter(const Vec3f &x) { mDetectorCenter = x; }

  float GetDetectorGain() const { return mDetectorGain; }
  float GetDetectorOffset() const { return mDetectorOffset; }
  void SetDetectorGain(float gain) { mDetectorGain = gain; }
  void SetDetectorOffset(float offset) { mDetectorOffset = offset; }

  Vec3f GetVolumeCenter() const { return mVolumeCenter; }


  // rotate a vector around z axis : by
  // [ [ cos -sin 0 ]
  //   [ sin  cos 0 ]
  //   [   0    0 1 ] ]
  template <typename T>
  static Vec3<T> RotateVector(const Vec3<T> &vec, float theta) {
    float c = cosf(theta);
    float s = sinf(theta);
    return Vec3<T>(c*vec[0] - s*vec[1],
                   s*vec[0] + c*vec[1],
                  vec[2]);
  }


private:

  // the location of the source in the rest position
  Vec3f mSourceLocation;
  Vec3f mSourceDirection;
  float mSourceCutoff;
  float mSourceFracAtCutoff;

  // the location of the detector in the rest position
  Vec3f mDetectorCenter;
  
  // detector orientation in the rest position
  Vec3f mDetectorForward;
  Vec3f mDetectorRight;
  Vec3f mDetectorUp;

  float mDetectorOffset;
  float mDetectorGain;
  
  // size of the detector in meters
  float mDetectorSize[2];

  // number of samples in the detector
  int mDetectorSamples[2];

  // size of each detector pixel in meters
  float mDetectorSampleSize[2];

  // length of the collimator - combined with pixel size gives collimation ratio
  float mDetectorCollimatorLength;
  float mDetectorCollimatorRatio; // based on pixel width, if height!=width

  // focal length of the collimator - positive values have focus behind the detector
  float mDetectorCollimatorFocalLength;

  // info about the projectections
  vector<float> mProjectionAngles;
  int mProjectionStride[3];

  // info about the volume we're trying the reconstruct
  Vec3f mVolumeCenter;
  Vec3f mVolumeSize;
  Vec3f mVolumeBounds[2];
  int mVolumeSamples[3];
  int mVolumeStride[3];
  int mVolumeNodeStride[3];
  
};

#endif
