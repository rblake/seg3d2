#ifndef _RBFInterface_H_
#define _RBFInterface_H_

#include <vector>
#include <string>

#include <rbf/Application/Tools/src/ScatteredData.h>
#include <rbf/Application/Tools/src/vec3.h>
#include <rbf/Application/Tools/src/RBF.h>
#include <rbf/Application/Tools/src/Surface.h>
//#include"SampleData.h"
//#include "horizon.h"

//using std::string; 
typedef vector< vector< vector<double> > > DataStructure;

class RBFInterface
{
public:
	ScatteredData *mySurfaceData;
	Surface *mySurface;
	RBF *mySurfaceRBF;
	Kernel myKernel;
	DataStructure value;

  // test
  size_t nx, ny, nz;
  double spacing_x, spacing_y, spacing_z;
  // test

//	void CreateSurface(string filename, vec3 myOrigin, vec3 mySize, vec3 mySampling);
	void CreateSurface(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing/*, vec3 mySampling*/);
//	RBFInterface(string filename, string dimensions);
	RBFInterface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing/*, vec3 mySampling*/);
//	RBFInterface();

  double getThresholdValue() const { return thresholdValue; }

private:
	void augmentNormalData(ScatteredData *data);
	vec3 findNormal(ScatteredData *data, int n);

  double thresholdValue;
};

#endif //_RBFInterface_H_ 
