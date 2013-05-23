
#ifndef MATERIAL_H
#define MATERIAL_H


#include <fstream>
#include <vector>
using std::vector;

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#define MATERIAL_ANGULAR_SAMPLES 1024


class Material {
public:

  bool InitFromFile(const char *fname, float lowKeV, float highKeV);
  bool InitFromString(const char *s, float lowKeV, float highKeV);

  float GetDensity() const { return mDensity; }
  float GetMassAttenuationCoefficient() const { return mMassAttenuationCoefficient; }
  float GetScatterFactor(float ncostheta) const;  // input 1-cos(theta)

  void GetScatterFactorArray(vector<float> &factors) const { factors = mScatterFactor; }

private:

  float ComputeCrossSection(const vector<float> &chi,
                            const vector<float> &F,
                            const vector<float> &S,
                            float theta,
                            float lowE, float highE) const;

  float Interpolate(const vector<float> &xs,
                    const vector<float> &ys,
                    float x) const;


  float mDensity;
  float mMassAttenuationCoefficient;
  float mMolecularWeight;
  vector<float> mScatterCrossSection;
  vector<float> mScatterFactor;
};


#endif
