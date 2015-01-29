//#include <string>
//
//#include "ScatteredData.h"
//#include "vec3.h"
//#include "SampleData.h"
//#include "RBF.h"
//#include "horizon.h"
//#include "fileIO.h"
#include <rbf/Application/Tools/src/RBFInterface.h>
#include <cmath>

// test
#include <fstream>
// test

//using std::string;

//void RBFInterface::CreateSurface(string filename, vec3 myOrigin, vec3 mySize, vec3 mySampling)
//{
//	mySurfaceData = new ScatteredData();
//	readSurfaceDataFile(filename, mySurfaceData);
//	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
//	//mySurfaceRBF->setDataReduction(Random);
//	myKernel = ThinPlate;
//	mySurface = new Surface(mySurfaceData, mySurfaceRBF);
//
//	//Construct RBFs
//	mySurface->computeRBF();
//
//	vec3 mySpacing(mySize[0]/mySampling[0], mySize[1]/mySampling[1], mySize[2]/mySampling[2]);
//	//printf("SPACING: %lf %lf %lf\n",mySpacing[0], mySpacing[1], mySpacing[2]);
//	value.resize((int)(mySampling[0]));
//	for(int i=0; i<mySampling[0]; i++)
//	{
//		//if(i%10==0)
//			//printf("%d/100 done\n", i); fflush(stdout);
//		value[i].resize((int)(mySampling[1]));
//		for(int j=0; j<mySampling[1]; j++)
//		{
//			//if(j%10==0)
//			//	printf("\t%d/100 done\n", j); fflush(stdout);
//			value[i][j].resize((int)(mySampling[2]));
//			for(int k=0; k<mySampling[2]; k++)
//			{
//				//if(k%10==0)
//				//	printf("\t\t%d/100 done\n", k); fflush(stdout);
//				vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
//				//std::cout<<"Computing Val ... "<<std::endl;
//				double myVal = mySurface->computeValue(location);
//				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
//				value[i][j][k]=myVal;
//			}
//		}
//	}
//}

//RBFInterface::RBFInterface(string filename, string dimensions)
//{
//	//read from the dimentsion file here	TODO
//	vec3 myOrigin(-30, -50, 80);
//	vec3 mySize(60, 50, 10);
//	vec3 mySampling(100, 100, 100);
//	CreateSurface(filename, myOrigin, mySize, mySampling);
//}
//
//RBFInterface::RBFInterface()
//{
//}

RBFInterface::RBFInterface(std::vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing/*, vec3 mySampling*/) : thresholdValue(0)
{
	CreateSurface(myData, myOrigin, mySize, mySpacing/*, mySampling*/);
}


void RBFInterface::CreateSurface(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySpacing/*, vec3 mySampling*/)
{
	vector<double> a,b,c,d;
	for(int i=0; i<myData.size(); i++)
	{
    std::cerr << "point : " << myData[i][0] << ", " << myData[i][1] << ", " << myData[i][2] << std::endl;
		a.push_back(myData[i][0]);
		b.push_back(myData[i][1]);
		c.push_back(myData[i][2]);
    // TODO: is this the value to threshold on?
    // If so, needs to be exposed in the interface!!!
		d.push_back(thresholdValue);
	}
  //This code figures out the bounds for the data selected

  std::vector<double>::iterator minx = std::min_element(a.begin(), a.end());
  std::vector<double>::iterator miny = std::min_element(b.begin(), b.end());
  std::vector<double>::iterator minz = std::min_element(c.begin(), c.end());
  std::vector<double>::iterator maxx = std::max_element(a.begin(), a.end());
  std::vector<double>::iterator maxy = std::max_element(b.begin(), b.end());
  std::vector<double>::iterator maxz = std::max_element(c.begin(), c.end());
  vec3 myMin(*minx, *miny, *minz), myMax(*maxx, *maxy, *maxz);
  myMin = myMin - 0.05*mySize;
  myMax = myMax + 0.05*mySize;

	mySurfaceData = new ScatteredData(a,b,c,d);
	augmentNormalData(mySurfaceData);
	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	mySurfaceRBF->setDataReduction(All);
  // TODO: let caller pick the kernel
	myKernel = ThinPlate;
	mySurface = new Surface(mySurfaceData, mySurfaceRBF);

	//Construct RBFs
	mySurface->computeRBF();

  //sanity check
  for(int i=0; i<mySurfaceData->fnc.size(); i++)
  {
    vec3 myLocation(mySurfaceData->x[0][i], mySurfaceData->x[1][i], mySurfaceData->x[2][i]);
    double myVal = mySurface->computeValue(myLocation);
    double error  = fabs(myVal - mySurfaceData->fnc[i]);
    if (error>1e-3)
    {
      printf("%lf\n", error);
      fflush(stdout);
    }
  }

  // TODO: dims and spacing of final dataset???
//	vec3 mySpacing(mySize[0]/mySampling[0], mySize[1]/mySampling[1], mySize[2]/mySampling[2]);

  // test
//  nx = mySampling[0];
//  ny = mySampling[1];
//  nz = mySampling[2];
//  spacing_x = mySpacing[0];
//  spacing_y = mySpacing[1];
//  spacing_z = mySpacing[2];
  // test

  // TODO: what happens when dims are changed to match input data?
  // TODO: use array (flat data structure) instead (cut down on data copying)?

	//printf("SPACING: %lf %lf %lf\n",mySpacing[0], mySpacing[1], mySpacing[2]);

  //Fill the values into the vector. In the first loop, we initialize the matrix with all values set to -100 (or - inf). In the second loop, we change the values from -100 to the correct value if the point in the domain described above.

  value.resize(static_cast<int>(mySize[0]));
  for(int i=0; i<mySize[0]; i++)
  {
    value[i].resize(static_cast<int>(mySize[1]));
    for(int j=0; j<mySize[1]; j++)
    {
      value[i][j].resize(static_cast<int>(mySize[2]), -100);
    }
  }

  for(int i=0; i<mySize[0]; i++)
  {
    vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX;
    if (location[0]<myMin[0]||location[0]>myMax[0])
      continue;
    for(int j=0; j<mySize[1]; j++)
    {
      location = myOrigin + mySpacing[1]*j*vec3::unitY;
      if (location[1]<myMin[1]||location[1]>myMax[1])
        continue;
      for(int k=0; k<mySize[2]; k++)
      {
        location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
        if (location[2]<myMin[2]||location[2]>myMax[2])
          continue;
        //std::cout<<"Computing Val ... "<<std::endl;
        double myVal = mySurface->computeValue(location);
        //printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
        value[i][j][k]=myVal;
      }
    }
  }

//  std::string filename = "surface.nrrd";
//  std::cout << "Writing file '" << filename << "'" << std::endl;
//  std::ofstream nrrd_file(filename.c_str(), std::ofstream::binary);
//
//  if (nrrd_file.is_open())
//  {
//    nrrd_file << "NRRD0001" << std::endl;
//    nrrd_file << "# Complete NRRD file format specification at:" << std::endl;
//    nrrd_file << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
//    nrrd_file << "type: float" << std::endl;
//    nrrd_file << "dimension: 3" << std::endl;
//    nrrd_file << "sizes: " << nx << " " << ny << " " << nz << std::endl;
//    nrrd_file << "axis mins: " << myOrigin[0] << ", " << myOrigin[1] << ", " << myOrigin[2] << std::endl;
//    nrrd_file << "spacings: " << spacing_x << " " << spacing_y << " " << spacing_z << std::endl;
//    nrrd_file << "centerings: cell cell cell" << std::endl;
//    nrrd_file << "endian: little" << std::endl;
//    nrrd_file << "encoding: raw" << std::endl;
//    nrrd_file << std::endl;
//
//    // write data portion
//    for(int k=0; k < nz; k++)
//    {
//      for(int j=0; j < ny; j++)
//      {
//        for(int i=0; i < nx; i++)
//        {
//          float val = value[i][j][k];
//          nrrd_file.write((char*)&val, sizeof(float));
//        }
//      }
//    }
//    nrrd_file.close();
//  }
}


vec3 RBFInterface::findNormal(ScatteredData *data, int n)
{
	int tot = data->x[0].size();
	int prev = (n-1)>=0?n-1:tot-1;
	int next = (n+1)<tot?n+1:0;

	while(data->x[2][prev]!=data->x[2][n])
	{
		prev = (prev-1)>=0?prev-1:tot-1;
	}

	while(data->x[2][next]!=data->x[2][n])
	{
		next = (next+1)<tot?next+1:0;
	}
	//printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);

	vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
	vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
	vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
	vec3 one = b-a;
	vec3 two = c-a;
	vec3 ret = one+two;
	return ret;
}

void RBFInterface::augmentNormalData(ScatteredData *data)
{
	int n = data->x[0].size();
	for(int i=0; i<n; i++)
	{
		vec3 myNormal = findNormal(data, i);
		myNormal = normalize(myNormal);
		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] + myNormal[j]);
		}
		data->fnc.push_back(10);

		for(int j=0; j<3; j++)
		{
			data->x[j].push_back(data->x[j][i] - myNormal[j]);
		}
		data->fnc.push_back(-10);
	}
}

