
#define _USE_MATH_DEFINES
#include <cmath>
#include <string.h>
#include <stdlib.h>

#ifndef WIN32
#define sscanf_s sscanf
#endif


#include "geometry.h"


bool Geometry::LoadFromFile(const char *fname) {
  std::ifstream f;
  f.open(fname);

  if (!f.is_open())
    return false;

  mProjectionAngles.clear();

  bool hasSourceLocation = false;
  bool hasSourceDirection = false;
  bool hasDetectorCenter = false;
  bool hasDetectorForward = false;
  bool hasDetectorRight = false;
  bool hasDetectorUp = false;
  bool hasDetectorSize = false;
  bool hasDetectorSamples = false;
  bool hasDetectorGain = false;
  bool hasDetectorOffset = false;
  bool hasDetectorCollimatorRatio = false;
  bool hasDetectorCollimatorFocalLength = false;
  bool hasProjectionAngles = false;
  bool hasVolumeSize = false;
  bool hasVolumeSamples = false;
  bool hasVolumeCenter = false;
  bool hasGPUBitfield = false;

  char line[1024];
  while (f.getline(line, 1024)) {
    if (strlen(line) < 1) continue;
    if (line[0] == '#') continue;

    if (!strchr(line, ':'))
      continue;

    char *key = line;
    char *value = strchr(line, ':');
    value[0] = 0;
    value++;

    if (strstr(key, "source location") &&
        sscanf_s(value, "%g, %g, %g", &mSourceLocation[0], &mSourceLocation[1], &mSourceLocation[2]) == 3) {
      hasSourceLocation = true;
    }

    if (strstr(key, "source direction") &&
        sscanf_s(value, "%g, %g, %g", &mSourceDirection[0], &mSourceDirection[1], &mSourceDirection[2]) == 3) {
      mSourceDirection.Normalize();
      hasSourceDirection = true;
    }

    if (strstr(key, "detector center") &&
        sscanf_s(value, "%g, %g, %g", &mDetectorCenter[0], &mDetectorCenter[1], &mDetectorCenter[2]) == 3) {
      hasDetectorCenter = true;
    }

    if (strstr(key, "detector forward") &&
        sscanf_s(value, "%g, %g, %g", &mDetectorForward[0], &mDetectorForward[1], &mDetectorForward[2]) == 3) {
      hasDetectorForward = true;
    }

    if (strstr(key, "detector right") &&
        sscanf_s(value, "%g, %g, %g", &mDetectorRight[0], &mDetectorRight[1], &mDetectorRight[2]) == 3) {
      hasDetectorRight = true;
    }

    if (strstr(key, "detector up") &&
        sscanf_s(value, "%g, %g, %g", &mDetectorUp[0], &mDetectorUp[1], &mDetectorUp[2]) == 3) {
      hasDetectorUp = true;
    }

    if (strstr(key, "detector size") &&
        sscanf_s(value, "%g, %g", &mDetectorSize[0], &mDetectorSize[1]) == 2) {
      hasDetectorSize = true;
    }

    if (strstr(key, "detector samples") &&
        sscanf_s(value, "%d, %d", &mDetectorSamples[0], &mDetectorSamples[1]) == 2) {
      hasDetectorSamples = true;
    }

    if (strstr(key, "detector gain") &&
        sscanf_s(value, "%g", &mDetectorGain) == 1) {
      hasDetectorGain = true;
    }

    if (strstr(key, "detector offset") &&
        sscanf_s(value, "%g", &mDetectorOffset) == 1) {
      hasDetectorOffset = true;
    }

    if (strstr(key, "detector collimator ratio") &&
        sscanf_s(value, "%g", &mDetectorCollimatorRatio) == 1) {
      hasDetectorCollimatorRatio = true;
    }

    if (strstr(key, "detector collimator focal length") &&
        sscanf_s(value, "%g", &mDetectorCollimatorFocalLength) == 1) {
      hasDetectorCollimatorFocalLength = true;
    }

    float pa;
    if (strstr(key, "projection angle") &&
        sscanf_s(value, "%g", &pa) == 1) {
      mProjectionAngles.push_back((float)(pa * M_PI / 180));
      hasProjectionAngles = true;
    }

    if (strstr(key, "volume size") &&
        sscanf_s(value, "%g, %g, %g", &mVolumeSize[0], &mVolumeSize[1], &mVolumeSize[2]) == 3) {
      hasVolumeSize = true;
    }

    if (strstr(key, "volume samples") &&
        sscanf_s(value, "%d, %d, %d", &mVolumeSamples[0], &mVolumeSamples[1], &mVolumeSamples[2]) == 3) {
      hasVolumeSamples = true;
    }

    if (strstr(key, "volume center") &&
        sscanf_s(value, "%g, %g, %g", &mVolumeCenter[0], &mVolumeCenter[1], &mVolumeCenter[2]) == 3) {
      hasVolumeCenter = true;
    }

    if (strstr(key, "use gpus") &&
        sscanf_s(value, "%d", &mGPUBitfield) == 1) {
      hasGPUBitfield = true;
    }
  }

  f.close();

  if (!hasSourceLocation)
    printf("error reading 'source location' from file\n");
  if (!hasSourceDirection)
    printf("error reading 'source direction' from file\n");
  if (!hasDetectorCenter)
    printf("error reading 'detector center' from file\n");
  if (!hasDetectorForward)
    printf("error reading 'detector forward' from file\n");
  if (!hasDetectorRight)
    printf("error reading 'detector right' from file\n");
  if (!hasDetectorUp)
    printf("error reading 'detector up' from file\n");
  if (!hasDetectorSize)
    printf("error reading 'detector size' from file\n");
  if (!hasDetectorSamples)
    printf("error reading 'detector samples' from file\n");
  if (!hasDetectorGain)
    printf("error reading 'detector gain' from file\n");
  if (!hasDetectorOffset)
    printf("error reading 'detector offset' from file\n");
  if (!hasDetectorCollimatorRatio)
    printf("error reading 'detector collimator ratio' from file\n");
  if (!hasDetectorCollimatorFocalLength)
    printf("error reading 'detector collimator focal length' from file\n");
  if (!hasProjectionAngles)
    printf("error reading 'projection angle' from file\n");
  if (!hasVolumeSize)
    printf("error reading 'volume size' from file\n");
  if (!hasVolumeSamples)
    printf("error reading 'volume samples' from file\n");
  if (!hasVolumeCenter)
    printf("error reading 'volume center' from file\n");
  if (!hasGPUBitfield)
    printf("error reading 'use gpus' from file\n");


  if (!(hasSourceLocation &&
        hasSourceDirection &&
        hasDetectorCenter &&
        hasDetectorForward &&
        hasDetectorRight &&
        hasDetectorUp &&
        hasDetectorSize &&
        hasDetectorSamples &&
        hasDetectorGain &&
        hasDetectorOffset &&
        hasDetectorCollimatorRatio &&
        hasDetectorCollimatorFocalLength &&
        hasProjectionAngles &&
        hasVolumeSize &&
        hasVolumeSamples &&
        hasVolumeCenter &&
        hasGPUBitfield))
    return false;


  mDetectorForward.Normalize();
  mDetectorRight.Normalize();
  mDetectorUp.Normalize();
  mDetectorSampleSize[0] = mDetectorSize[0] / mDetectorSamples[0];
  mDetectorSampleSize[1] = mDetectorSize[1] / mDetectorSamples[1];

  mProjectionStride[0] = 1;
  mProjectionStride[1] = mProjectionStride[0] * mDetectorSamples[0];
  mProjectionStride[2] = mProjectionStride[1] * mDetectorSamples[1];

  mVolumeStride[0] = 1;
  mVolumeStride[1] = mVolumeStride[0] * mVolumeSamples[0];
  mVolumeStride[2] = mVolumeStride[1] * mVolumeSamples[1];

  mVolumeNodeStride[0] = 1;
  mVolumeNodeStride[1] = mVolumeNodeStride[0] * (mVolumeSamples[0]+1);
  mVolumeNodeStride[2] = mVolumeNodeStride[1] * (mVolumeSamples[1]+1);

  mVolumeBounds[0] = mVolumeCenter - 0.5 * mVolumeSize;
  mVolumeBounds[1] = mVolumeCenter + 0.5 * mVolumeSize;

  mDetectorCollimatorLength = mDetectorSampleSize[0] * mDetectorCollimatorRatio;

  WorldToVolume(mSourceLocation, mSourceLocationVol);


  return true;
}


bool Geometry::SaveToFile(const char *fname) const {
  std::ofstream f;
  f.open(fname);

  if (!f.is_open())
    return false;

  f << "source location: " << mSourceLocation[0] << ", " << mSourceLocation[1] << ", " << mSourceLocation[2] << std::endl;
  f << "source direction: " << mSourceDirection[0] << ", " << mSourceDirection[1] << ", " << mSourceDirection[2] << std::endl;
  f << "detector center: " << mDetectorCenter[0] << ", " << mDetectorCenter[1] << ", " << mDetectorCenter[2] << std::endl;
  f << "detector forward: " << mDetectorForward[0] << ", " << mDetectorForward[1] << ", " << mDetectorForward[2] << std::endl;
  f << "detector right: " << mDetectorRight[0] << ", " << mDetectorRight[1] << ", " << mDetectorRight[2] << std::endl;
  f << "detector up: " << mDetectorUp[0] << ", " << mDetectorUp[1] << ", " << mDetectorUp[2] << std::endl;
  f << "detector size: " << mDetectorSize[0] << ", " << mDetectorSize[1] << std::endl;
  f << "detector samples: " << mDetectorSamples[0] << ", " << mDetectorSamples[1] << std::endl;
  f << "detector gain: " << mDetectorGain << std::endl;
  f << "detector offset: " << mDetectorOffset << std::endl;
  f << "detector collimator ratio: " << mDetectorCollimatorRatio << std::endl;
  f << "detector collimator focal length: " << mDetectorCollimatorFocalLength << std::endl;
  f << "volume size: " << mVolumeSize[0] << ", " << mVolumeSize[1] << ", " << mVolumeSize[2] << std::endl;
  f << "volume samples: " << mVolumeSamples[0] << ", " << mVolumeSamples[1] << ", " << mVolumeSamples[2] << std::endl;
  f << "volume center: " << mVolumeCenter[0] << ", " << mVolumeCenter[1] << ", " << mVolumeCenter[2] << std::endl;

  f << "use gpus: " << mGPUBitfield << std::endl;

  for (size_t i=0; i<mProjectionAngles.size(); i++) {
    f << "projection angle: " << (mProjectionAngles[i]*180/M_PI) << std::endl;
  }


  f.close();

  return true;
}


// convert projection index types
int Geometry::ProjectionCoordToIndex(int p, int x, int y) const {
  return mProjectionStride[0]*x + mProjectionStride[1]*y + mProjectionStride[2]*p;
}
void Geometry::ProjectionIndexToCoord(int i, int &p, int &x, int &y) const {
  p = i/mProjectionStride[2];
  i %= mProjectionStride[2];

  y = i/mProjectionStride[1];
  i %= mProjectionStride[1];

  x = i/mProjectionStride[0];
}


// convert volume index types
int Geometry::VolumeCoordToIndex(int x, int y, int z) const {
  return mVolumeStride[0]*x + mVolumeStride[1]*y + mVolumeStride[2]*z;
}
void Geometry::VolumeIndexToCoord(int i, int &x, int &y, int &z) const {
  z = i/mVolumeStride[2];
  i %= mVolumeStride[2];

  y = i/mVolumeStride[1];
  i %= mVolumeStride[1];

  x = i/mVolumeStride[0];
}


int Geometry::VolumeNodeCoordToIndex(int x, int y, int z) const {
  return mVolumeNodeStride[0]*x + mVolumeNodeStride[1]*y + mVolumeNodeStride[2]*z;
}
void Geometry::VolumeIndexToNodeCoord(int i, int &x, int &y, int &z) const {
  z = i/mVolumeNodeStride[2];
  i %= mVolumeNodeStride[2];

  y = i/mVolumeNodeStride[1];
  i %= mVolumeNodeStride[1];

  x = i/mVolumeNodeStride[0];
}


// map point from world coordinates to volume pixel coordinates, return false if outside the volume.
bool Geometry::WorldToVolume(const Vec3f &world, Vec3f &volume) const {
  bool inside = true;
  for (int i=0; i<3; i++) {
    volume[i] = mVolumeSamples[i] * (world[i]-mVolumeBounds[0][i]) / (mVolumeBounds[1][i]-mVolumeBounds[0][i]);
    if (volume[i] < 0 || volume[i] >= mVolumeSamples[i])
      inside = false; // outside the volume
  }
  return inside;
}

// map point from volume pixel coordinates to world coordinates
void Geometry::VolumeToWorld(const Vec3f &volume, Vec3f &world) const {
  for (int i=0; i<3; i++) {
    world[i] = mVolumeBounds[0][i] + (volume[i] / mVolumeSamples[i]) * (mVolumeBounds[1][i]-mVolumeBounds[0][i]);
  }
}

// map direction from world coordinates to volume pixel coordinates
void Geometry::WorldToVolumeDir(const Vec3f &world, Vec3f &volume) const {
  for (int i=0; i<3; i++) {
    volume[i] = world[i] * (mVolumeSamples[i] / (mVolumeBounds[1][i]-mVolumeBounds[0][i]));
  }
}


// return the min and max distance along the ray where it intersects with the volume's bounding box.
// tmin/tmax may be negative!  Return false if no intersection.
bool Geometry::RayVolumeIntersection(const Vec3f &origin, const Vec3f &dir,
                                     float &tmin, float &tmax) const {
  float tymin, tymax, tzmin, tzmax;
  Vec3f inv_direction(1/dir[0], 1/dir[1], 1/dir[2]);
  int sign[3] = { inv_direction[0]<0, inv_direction[1]<0, inv_direction[2]<0 };

  tmin = (mVolumeBounds[sign[0]][0] - origin[0]) * inv_direction[0];
  tmax = (mVolumeBounds[1-sign[0]][0] - origin[0]) * inv_direction[0];
  tymin = (mVolumeBounds[sign[1]][1] - origin[1]) * inv_direction[1];
  tymax = (mVolumeBounds[1-sign[1]][1] - origin[1]) * inv_direction[1];
  if ( (tmin > tymax) || (tymin > tmax) ) 
    return false;
  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;
  tzmin = (mVolumeBounds[sign[2]][2] - origin[2]) * inv_direction[2];
  tzmax = (mVolumeBounds[1-sign[2]][2] - origin[2]) * inv_direction[2];
  if ( (tmin > tzmax) || (tzmin > tmax) ) 
    return false;
  if (tzmin > tmin)
    tmin = tzmin;
  if (tzmax < tmax)
    tmax = tzmax;
  return true;
}


// for a point in world space, find the voxels and weights for trilinear sampling at that point.
// return false if outside of the volume
bool Geometry::WorldToTrilinearWeights(const Vec3f &world, int indices[8], float weights[8]) const {
  Vec3f volume;
  if (!WorldToVolume(world, volume))
    return false;

  int vi[3] = { (int)volume[0], (int)volume[1], (int)volume[2] };
  float frac[3] = { volume[0]-vi[0], volume[1]-vi[1], volume[2]-vi[2] };

  float xfracs[2] = { 1-frac[0], frac[0] };
  float yfracs[2] = { 1-frac[1], frac[1] };
  float zfracs[2] = { 1-frac[2], frac[2] };

  for (int x=0; x<2; x++) {
    for (int y=0; y<2; y++) {
      for (int z=0; z<2; z++) {

        int i = x*4 + y*2 + z;
        indices[i] = VolumeNodeCoordToIndex(vi[0]+x, vi[1]+y, vi[2]+z);
        weights[i] = xfracs[x]*yfracs[y]*zfracs[z];
      }
    }
  }

  return true;
}


Vec3f Geometry::GetSourcePosition() const {
  return mSourceLocation;
}

Vec3f Geometry::GetSourceDirection() const {
  return mSourceDirection;
}


void Geometry::SetSourceAttenMap(const float *sourceAttenMap,
                                 int width, int height) {
  mSourceAttenMapSize[0] = width;
  mSourceAttenMapSize[1] = height;

  mSourceAttenMap.resize(width*height);
  for (int i=0; i<width*height; i++)
    mSourceAttenMap[i] = sourceAttenMap[i];
}


float Geometry::GetSourceIntensityThroughPoint(const Vec3f &pos, int proj) const {

  Vec3<float> ipos = GetInverseRotatedPoint(pos, proj);

  // just use nearest sample
  // projection transform asuming detector was 0.19m from source and is 0.205m on a side
  // x axis is swapped, but y is correct
  int sx = (int)((-pos[0] / (pos[2]-mSourceLocation[2])) * mSourceAttenMapSize[0] * (0.19 / 0.205) + mSourceAttenMapSize[0]*0.5 + 0.5f);
  int sy = (int)(( pos[1] / (pos[2]-mSourceLocation[2])) * mSourceAttenMapSize[1] * (0.19 / 0.205) + mSourceAttenMapSize[1]*0.5 + 0.5f);

  if (sx<0) sx = 0;
  else if (sx>=mSourceAttenMapSize[0]) sx = mSourceAttenMapSize[0]-1;

  if (sy<0) sy = 0;
  else if (sy>=mSourceAttenMapSize[1]) sy = mSourceAttenMapSize[1]-1;

  return mSourceAttenMap[sy*mSourceAttenMapSize[0] + sx];
}


void Geometry::SetVoxelSize(const Vec3f &size) {
  for (int i=0; i<3; i++) {
    mVolumeSamples[i] = (int)(((mVolumeBounds[1][i] - mVolumeBounds[0][i]) / size[i]) + 0.5);
  }

  mVolumeStride[0] = 1;
  mVolumeStride[1] = mVolumeStride[0] * mVolumeSamples[0];
  mVolumeStride[2] = mVolumeStride[1] * mVolumeSamples[1];

  mVolumeNodeStride[0] = 1;
  mVolumeNodeStride[1] = mVolumeNodeStride[0] * (mVolumeSamples[0]+1);
  mVolumeNodeStride[2] = mVolumeNodeStride[1] * (mVolumeSamples[1]+1);
}


void Geometry::GetDetectorRay(int p, int x, int y, bool centered,
                              Vec3f &origin, Vec3f &dir) const {

  Vec3f pixelSample, colSample;

  float focalScale;
  if (mDetectorCollimatorFocalLength == 0)
    focalScale = 1;
  else
    focalScale = (mDetectorCollimatorFocalLength + mDetectorCollimatorLength) / mDetectorCollimatorFocalLength;

  if (centered) {
    pixelSample = (mDetectorCenter + 
                   mDetectorRight * (mDetectorSampleSize[0] * (x - mDetectorSamples[0]*0.5f + 0.5f)) +
                   mDetectorUp * (mDetectorSampleSize[1] * (y - mDetectorSamples[1]*0.5f + 0.5f)));
    colSample = (mDetectorCenter + 
                 mDetectorRight * (focalScale * mDetectorSampleSize[0] * (x - mDetectorSamples[0]*0.5f + 0.5f)) +
                 mDetectorUp * (focalScale * mDetectorSampleSize[1] * (y - mDetectorSamples[1]*0.5f + 0.5f)) +
                 mDetectorForward * mDetectorCollimatorLength);
  }
  else {
    pixelSample = (mDetectorCenter + 
                   mDetectorRight * (mDetectorSampleSize[0] * (x - mDetectorSamples[0]*0.5f + ((float)rand() / RAND_MAX))) +
                   mDetectorUp * (mDetectorSampleSize[1] * (y - mDetectorSamples[1]*0.5f + ((float)rand() / RAND_MAX))));
    colSample = (mDetectorCenter + 
                 mDetectorRight * (focalScale * mDetectorSampleSize[0] * (x - mDetectorSamples[0]*0.5f + ((float)rand() / RAND_MAX))) +
                 mDetectorUp * (focalScale * mDetectorSampleSize[1] * (y - mDetectorSamples[1]*0.5f + ((float)rand() / RAND_MAX))) +
                 mDetectorForward * mDetectorCollimatorLength);
  }

  origin = RotateVector(pixelSample, mProjectionAngles[p]);
  dir = RotateVector((colSample-pixelSample).Normalized(), mProjectionAngles[p]);
}


void Geometry::ProjectionPixelToPhysicalOffset(float px, float py, float &mx, float &my) const {
  mx = (px - mDetectorSamples[0]*0.5f) * mDetectorSampleSize[0];
  my = (py - mDetectorSamples[1]*0.5f) * mDetectorSampleSize[1];
}

void Geometry::ProjectionPhysicalOffsetToPixel(float mx, float my, float &px, float &py) const {
  px = mx/mDetectorSampleSize[0] + mDetectorSamples[0]*0.5f;
  py = my/mDetectorSampleSize[1] + mDetectorSamples[1]*0.5f;
}
