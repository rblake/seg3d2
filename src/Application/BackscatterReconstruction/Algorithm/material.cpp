
#include <string.h>

#include "material.h"

#ifndef WIN32
#define sscanf_s sscanf
#endif


template <typename T>
T square(T x) { return x*x; }


bool Material::InitFromFile(const char *fname, float lowKeV, float highKeV) {

  std::ifstream file;
  file.open(fname);

  if (!file.is_open()) {
    printf("ScatterDistribution: error opening file %s\n", fname);
    return false;
  }

  bool hasMolecularWeight = false;
  bool hasScatterFactors = false;
  bool hasAttenuation = false;
  bool hasDensity = false;
  

  vector<float> chi, F, S;

  char line[1024];
  while (file.getline(line, 1024)) {
    if (strlen(line) < 1) continue;
    if (line[0] == '#') continue;

    if (!strchr(line, ':'))
      continue;

    char *key = line;
    char *value = strchr(line, ':');
    value[0] = 0;
    value++;

    float tchi, tF, tS;
    if (strstr(key, "scatter") &&
        sscanf_s(value, "%g, %g, %g", &tchi, &tF, &tS) == 3) {
      chi.push_back(tchi);
      F.push_back(tF);
      S.push_back(tS);
      hasScatterFactors = true;
    }

    if (strstr(key, "molecular weight") &&
        sscanf_s(value, "%g", &mMolecularWeight) == 1) {
      hasMolecularWeight = true;
    }

    if (strstr(key, "mass attenuation coefficient") &&
        sscanf_s(value, "%g", &mMassAttenuationCoefficient) == 1) {
      hasAttenuation = true;
    }

    if (strstr(key, "density") &&
        sscanf_s(value, "%g", &mDensity) == 1) {
      hasDensity = true;
    }
  }
  file.close();


  // make sure we read everything
  if (!hasScatterFactors)
    printf("error reading 'scatter' from material file\n");
  if (!hasMolecularWeight)
    printf("error reading 'molecular weight' from material file\n");
  if (!hasAttenuation)
    printf("error reading 'mass attenuation coefficient' from material file\n");
  if (!hasDensity)
    printf("error reading 'density' from material file\n");

  if (!(hasScatterFactors &&
        hasMolecularWeight &&
        hasAttenuation &&
        hasDensity))
    return false;


  // resample to have a uniformly sampled scatter cross section
  int numAngularSamples = MATERIAL_ANGULAR_SAMPLES;
  mScatterCrossSection.clear();
  mScatterFactor.clear();

  for (int thetai=0; thetai<numAngularSamples; thetai++) {
    float theta = (float)M_PI * thetai/(numAngularSamples-1.0f);
    mScatterCrossSection.push_back(ComputeCrossSection(chi, F, S, theta, lowKeV, highKeV));
  }

  for (int i=0; i<numAngularSamples; i++) {
    float ncostheta = 1 - 2*(float)i/(numAngularSamples-1);
    float theta = acos(ncostheta);
    float avagadro = 6.02214179e23f;
    if (mMolecularWeight > 0)
      mScatterFactor.push_back(ComputeCrossSection(chi, F, S, theta, lowKeV, highKeV) *
                               (avagadro / mMolecularWeight) * mDensity);
    else
      mScatterFactor.push_back(0);
  }

  return true;
}


float Material::GetScatterFactor(float ncostheta) const {
  float cthetas = (ncostheta * (mScatterCrossSection.size()-1));
  int cthetai = (int)cthetas;
  float cthetaf = cthetas-cthetai;

  if (cthetai<0)
    return mScatterFactor.front();
  else if (cthetai>=mScatterFactor.size()-1)
    return mScatterFactor.back();
  else
    return (1-cthetaf)*mScatterFactor[cthetai] + cthetaf*mScatterFactor[cthetai+1];
}


float Material::ComputeCrossSection(const vector<float> &chi,
                                    const vector<float> &F,
                                    const vector<float> &S,
                                    float theta,
                                    float lowKeV, float highKeV) const {

  float costheta = cosf(theta);

  float averageCrossSection = 0;
  int numEnergySamples = 100;
  for (int ei=0; ei<numEnergySamples; ei++) {
    float energy = lowKeV + ei/(numEnergySamples-1.0f) * (highKeV-lowKeV);
    float lambda = 12.39852f / energy;
    float k = energy / 511.0034f;

    float sampleChi = sin(theta/2) / lambda;
    float sampleF = Interpolate(chi, F, sampleChi);
    float sampleS = Interpolate(chi, S, sampleChi);

    float re = 2.817938e-15f;
    float coherent_dsigma_domega = (re*re/2) * (1+square(costheta));
    float incoherent_dsigma_domega = ((re*re/2) * 
                                      (1/square(1+k*(1-costheta))) *
                                      (1+square(costheta)+square(k)*square(1-costheta)) / (1+k*(1-costheta)));


    float coherent_cross = sampleF*sampleF*coherent_dsigma_domega;
    float incoherent_cross = sampleS*incoherent_dsigma_domega;
    averageCrossSection += coherent_cross + incoherent_cross;
  }

  return (averageCrossSection / numEnergySamples);
}


float Material::Interpolate(const vector<float> &xs,
                            const vector<float> &ys,
                            float x) const {


  int lowI = 0;
  int highI = 0;

  // extrapolate to low values
  if (x <= xs.front()) {
    lowI = 0;
    highI = 1;
  }

  // extrapolate to high values
  else if (x >= xs.back()) {
    for (highI=1; highI<(int)xs.size(); highI++) {
      if (xs[highI] > x) {
        lowI = highI-1;
        break;
      }
    }
  }

  // linear interpolate
  else {
    lowI = (int)xs.size()-2;
    highI = lowI+1;
  }


  float alpha = (x-xs[lowI]) / (xs[highI]-xs[lowI]);
  return ys[lowI] + alpha*(ys[highI]-ys[lowI]);
}

