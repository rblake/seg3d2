#include <string>

#include "ScatteredData.h"
#include "vec3.h"
#include "SampleData.h"
#include "RBF.h"
#include "horizon.h"
#include "fileIO.h"
#include "API.h"

using std::string;

void API::CreateSurface(string filename, vec3 myOrigin, vec3 mySize, vec3 mySampling)
{
	mySurfaceData = new ScatteredData();
	readSurfaceDataFile(filename, mySurfaceData);
	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	//mySurfaceRBF->setDataReduction(Random);
	myKernel = ThinPlate;
	mySurface = new Horizon(mySurfaceData, mySurfaceRBF);

	//Construct RBFs
	mySurface->computeRBF();

	vec3 mySpacing(mySize[0]/mySampling[0], mySize[1]/mySampling[1], mySize[2]/mySampling[2]);
	//printf("SPACING: %lf %lf %lf\n",mySpacing[0], mySpacing[1], mySpacing[2]);
	value.resize((int)(mySampling[0]));
	for(int i=0; i<mySampling[0]; i++)
	{
		//if(i%10==0)
			//printf("%d/100 done\n", i); fflush(stdout);
		value[i].resize((int)(mySampling[1]));
		for(int j=0; j<mySampling[1]; j++)
		{
			//if(j%10==0)
			//	printf("\t%d/100 done\n", j); fflush(stdout);
			value[i][j].resize((int)(mySampling[2]));
			for(int k=0; k<mySampling[2]; k++)
			{
				//if(k%10==0)
				//	printf("\t\t%d/100 done\n", k); fflush(stdout);
				vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
				//std::cout<<"Computing Val ... "<<std::endl;
				double myVal = mySurface->computeValue(location);
				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
				value[i][j][k]=myVal;
			}
		}
	}
}

API::API(string filename, string dimensions)
{
	//read from the dimentsion file here	TODO
	vec3 myOrigin(-30, -50, 80);
	vec3 mySize(60, 50, 10);
	vec3 mySampling(100, 100, 100);
	CreateSurface(filename, myOrigin, mySize, mySampling);
}

API::API()
{
}

API::API(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling)
{
	CreateSurface(myData, myOrigin, mySize, mySampling);
}

void API::CreateSurface(vector<vec3> myData, vec3 myOrigin, vec3 mySize, vec3 mySampling)
{
	vector<double> a,b,c,d;
	for(int i=0; i<myData.size(); i++)
	{
		a.push_back(myData[i][0]);
		b.push_back(myData[i][1]);
		c.push_back(myData[i][2]);
		d.push_back(0);
	}
	mySurfaceData = new ScatteredData(a,b,c,d);
	mySurfaceRBF = new RBF(mySurfaceData, myKernel);
	//mySurfaceRBF->setDataReduction(Random);
	myKernel = ThinPlate;
	mySurface = new Horizon(mySurfaceData, mySurfaceRBF);

	//Construct RBFs
	mySurface->computeRBF();

	vec3 mySpacing(mySize[0]/mySampling[0], mySize[1]/mySampling[1], mySize[2]/mySampling[2]);
	//printf("SPACING: %lf %lf %lf\n",mySpacing[0], mySpacing[1], mySpacing[2]);
	value.resize((int)(mySampling[0]));
	for(int i=0; i<mySampling[0]; i++)
	{
		//if(i%10==0)
			//printf("%d/100 done\n", i); fflush(stdout);
		value[i].resize((int)(mySampling[1]));
		for(int j=0; j<mySampling[1]; j++)
		{
			//if(j%10==0)
			//	printf("\t%d/100 done\n", j); fflush(stdout);
			value[i][j].resize((int)(mySampling[2]));
			for(int k=0; k<mySampling[2]; k++)
			{
				//if(k%10==0)
				//	printf("\t\t%d/100 done\n", k); fflush(stdout);
				vec3 location = myOrigin + mySpacing[0]*i*vec3::unitX + mySpacing[1]*j*vec3::unitY + mySpacing[2]*k*vec3::unitZ;
				//std::cout<<"Computing Val ... "<<std::endl;
				double myVal = mySurface->computeValue(location);
				//printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
				value[i][j][k]=myVal;
			}
		}
	}
}
