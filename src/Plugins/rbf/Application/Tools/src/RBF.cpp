//-------------------------------------------------------------------
//
//  Permission is  hereby  granted, free  of charge, to any person
//  obtaining a copy of this software and associated documentation
//  files  ( the "Software" ),  to  deal in  the  Software without
//  restriction, including  without limitation the rights to  use,
//  copy, modify,  merge, publish, distribute, sublicense,  and/or
//  sell copies of the Software, and to permit persons to whom the
//  Software is  furnished  to do  so,  subject  to  the following
//  conditions:
//
//  The above  copyright notice  and  this permission notice shall
//  be included  in  all copies  or  substantial  portions  of the
//  Software.
//
//  THE SOFTWARE IS  PROVIDED  "AS IS",  WITHOUT  WARRANTY  OF ANY
//  KIND,  EXPRESS OR IMPLIED, INCLUDING  BUT NOT  LIMITED  TO THE
//  WARRANTIES   OF  MERCHANTABILITY,  FITNESS  FOR  A  PARTICULAR
//  PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL THE AUTHORS OR
//  COPYRIGHT HOLDERS  BE  LIABLE FOR  ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//  USE OR OTHER DEALINGS IN THE SOFTWARE.
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#include <rbf/Application/Tools/src/ScatteredData.h>
#include <rbf/Application/Tools/src/RBF.h>
#include <rbf/Application/Tools/src/SparseMatrix.h>
#include <rbf/Application/Tools/src/vec3.h>
#include <rbf/Application/Tools/src/LinearSolver.h>
#include <rbf/Application/Tools/src/FMM.h>

#include <cmath>
#include <cstdio>

// STL Includes
#include <vector>
#include <utility>
#include <algorithm>


//RBF::RBF()
//{
//}

RBF::RBF(ScatteredData *myData, Kernel myKernel) :
  kernel(ThinPlate), acceleration(None), dataReduction(All)
{
	setKernel(myKernel);
	setData(myData);
}

void RBF::setKernel(Kernel myKernel)
{
	kernel = myKernel;
}

void RBF::setData(ScatteredData *myData)
{
	completeData = myData;
}

void RBF::setAcceleration(Acceleration myAcceleration)
{
	acceleration = myAcceleration;
	switch(acceleration)
	{
		case FastMultipole:
			fmm = new FMM();
			break;
	}
}
void RBF::setDataReduction(DataReduction myDataReduction)
{
	dataReduction = myDataReduction;
}

void RBF::computeFunction()
{
  const int ERROR_THRESHOLD = 20; // 20 or 2
	data = new ScatteredData;
	switch(dataReduction)
	{
		case All:
			data->setData(completeData->x[0], completeData->x[1], completeData->x[2], completeData->fnc);
			computeFunctionForData();
			break;
      
		case Random:
			bool *added;
			int n = completeData->fnc.size();
			printf("%d\n",n);
			added = new bool[n];
      
			for(int i=0; i<n; i++)
				added[i]=false;
      
			for(int i=0; i<25; i++)
			{
				int j = rand()%n;
				if(added[j])
				{
					i--;
					continue;
				}
				added[j]=true;
				data->x[0].push_back(completeData->x[0][j]);
				data->x[1].push_back(completeData->x[1][j]);
				data->x[2].push_back(completeData->x[2][j]);
				data->fnc.push_back(completeData->fnc[j]);
				printf("%d %lf %lf %lf %lf\n", j, completeData->x[0][j],completeData->x[1][j],completeData->x[2][j],completeData->fnc[j]);
			}
			std::vector<std::pair<double, int> > error;
			bool smallError=false;
			while(!smallError)
			{
				computeFunctionForData();
				computeErrorForData(error);
				std::sort(error.begin(), error.end());
				std::reverse(error.begin(), error.end());
				
				printf("Largest error: %lf\n", error[0].first);
        if(error[0].first < ERROR_THRESHOLD)
				{
					smallError=true;
					continue;
				}
				for(int k=0; k<5; k++)
				{
					//printf("Error %d: %lf\n", k, error[k].first);
					if(error[k].first>ERROR_THRESHOLD || error[k].first!=error[k].first)
					{
						int j = error[k].second;
						printf("Adding data point %d\n", j);
						added[j]=true;
						data->x[0].push_back(completeData->x[0][j]);
						data->x[1].push_back(completeData->x[1][j]);
						data->x[2].push_back(completeData->x[2][j]);
						data->fnc.push_back(completeData->fnc[j]);
					}
				}
        
			}
			printf("Total no. of data point: %d\n",  data->fnc.size()); fflush(stdout);
	}
}
void RBF::computeFunctionForData()
{
	switch(acceleration)
	{
		case FastMultipole:
			fmmComputeFunction();
		case None:
		default:
			int n = data->fnc.size();
			printf("Solving linear equations: \n"); fflush(stdout);
			LinearSolver rbfSolver;
			SparseMatrix rbfMatrix(n);
			printf("Constructing matrix ... "); fflush(stdout);
			for(int i=0; i<n; i++)
			{
				for(int j=0; j<n; j++)
				{
					//printf("%d %d ", i,j); fflush(stdout);
					double val = computeKernel(i,j);
					//printf("%lf\n", val); fflush(stdout);
					rbfMatrix.push_back(i,j,val);
				}
			}
			printf("Done\n"); fflush(stdout);
			rbfSolver.setMatrix(&rbfMatrix);
			printf("Running BiCGSTAB Iterations ... "); fflush(stdout);
			rbfSolver.biCGStab(data->fnc, coeff);
			printf("Done\n"); fflush(stdout);
	}
}

double RBF::computeValue(vec3 x)
{
	switch(acceleration)
	{
		case FastMultipole:
			return fmmComputeValue(x);
		case None:
		default:
      double sum=0;
      for(int i=0; i<coeff.size(); i++)
        sum+=coeff[i]*computeKernel(i, x);
      return sum;
	}
}

void RBF::computeErrorForData(std::vector<std::pair<double, int> > &error)
{
	int n = completeData->fnc.size();
	error.clear();
	for (int i=0; i<n; i++)
	{
		vec3 x(completeData->x[0][i],completeData->x[1][i],completeData->x[2][i]);
		double err = completeData->fnc[i]-computeValue(x);
		error.push_back(std::make_pair(err,i));
	}
}

double RBF::computeKernel(int i, int j)
{
	double r = sqrt( (data->x[0][i] - data->x[0][j])*(data->x[0][i] - data->x[0][j]) +
                  (data->x[1][i] - data->x[1][j])*(data->x[1][i] - data->x[1][j]) +
                  (data->x[2][i] - data->x[2][j])*(data->x[2][i] - data->x[2][j]) );
  
	return computeRadialFunction(r);
  
}

double RBF::computeKernel(int i, vec3 b)
{
	double r = sqrt( (data->x[0][i] - b[0])*(data->x[0][i] - b[0]) +
                  (data->x[1][i] - b[1])*(data->x[1][i] - b[1]) +
                  (data->x[2][i] - b[2])*(data->x[2][i] - b[2]) );
  
	return computeRadialFunction(r);
}

double RBF::computeRadialFunction(double r)
{
	double c = 0.1;
	switch(kernel)
	{
		case Guassian:
			r = r*0.01;
			return 1.0/sqrt(r*r + c*c);
			break;
		case ThinPlate:
			return r*r*log(r+c);
			break;
		case MultiQuadratic:
			return sqrt(r*r + c*c);
		default:
			return r;
			break;
	}
	return 0;
}


//FMM Codes
void RBF::fmmComputeFunction()
{
	fmmBuildTree();
}

void RBF::fmmBuildTree()
{
	std::vector<int> myIndices;
	int n = data->x[0].size();
	for(int i=0; i<n; i++)
		myIndices.push_back(i);
	fmm->tree = new BHNode();
	fmm->tree->box.min = vec3::zero;
	vec3 corner(100,100,100);
	fmm->tree->box.max = corner;
	fmmBuildTree(myIndices, fmm->tree);
	//printf("Tree Built\n");
	//fmmPrintTree(fmm->tree, 0);
}

void RBF::fmmPrintTree(BHNode *myNode, int stack)
{
	if (stack>4)
		return;
	for(int j=0; j<stack; j++)
	{
		printf(" ");
	}
	printf("%d %d %d\n", myNode, myNode->index, myNode->pts.size());
	for(int i=0; i<8; i++)
	{
		if(myNode->nodes[i]!=0)
			fmmPrintTree(myNode->nodes[i], stack+1);
	}
}

void RBF::fmmBuildTree(std::vector<int> &myPoints, BHNode *myNode)
{
  
	//printf("[%lf %lf %lf] [%lf %lf %lf] %d\n", myNode->box.min[0], myNode->box.min[1], myNode->box.min[2], myNode->box.max[0], myNode->box.max[1], myNode->box.max[2], myPoints.size());
	std::vector<int> children[8];
	int n = myPoints.size();
  
	myNode->index = fmm->numOfNodes;
	fmm->numOfNodes += 1;
	fmm->nodePointer.push_back(myNode);
  
	myNode->mass = n;
	myNode->leaf = (n<=1)?true:false;
	myNode->center = vec3::zero;
  
	myNode->coeff = 1; //REPLACE
	//add all the coefficients
  
  
	for(int i=0; i<myPoints.size(); i++)
		myNode->pts.push_back(myPoints[i]);
  
	for(int i=0; i<n; i++)
	{
		vec3 location(data->x[0][myPoints[i]], data->x[1][myPoints[i]],data->x[2][myPoints[i]]);
		myNode->center = myNode->center + (location/n);
	}
  
	if (n==1)
	{
		for(int i=0; i<8; i++)
			myNode->nodes[i]=0;
		return;
	}
	if(length(myNode->box.max - myNode->box.min) < 1e-6)
	{
		for(int i=0; i<8; i++)
			myNode->nodes[i]=0;
		return;
	}
  
	vec3 mid((myNode->box.getMin()+myNode->box.getMax())/2);
	for(int i=0; i<n; i++)
	{
		int octant=0;
		//FIND OCTANTS
		if(data->x[0][myPoints[i]] > mid[0])
			octant+=1;
		if(data->x[1][myPoints[i]] > mid[1])
			octant+=2;
		if(data->x[2][myPoints[i]] > mid[2])
			octant+=4;
    
		//printf("%d %d %d %lf %lf %lf\n", i,octant, myPoints[i], data->x[0][myPoints[i]],data->x[1][myPoints[i]], data->x[2][myPoints[i]]);
    
		children[octant].push_back(myPoints[i]);
	}
	for(int i=0; i<8; i++)
	{
		if (children[i].size() == 0)
		{
			myNode->nodes[i] = 0;
			continue;
		}
    
		myNode->nodes[i] = new BHNode();
    
		int hash = i;
		for(int j=0; j<3; j++)
		{
			if(hash%2==0)
			{
				myNode->nodes[i]->box.min[j]=myNode->box.min[j];
				myNode->nodes[i]->box.max[j]=mid[j];
			}
			else
			{
				myNode->nodes[i]->box.min[j]=mid[j];
				myNode->nodes[i]->box.max[j]=myNode->box.max[j];
			}
			hash /= 2;
		}
    
		
		fmmBuildTree(children[i], myNode->nodes[i]);
	}
}

double RBF::fmmComputeValue(vec3 x)
{
	double val = fmmComputeValueRecurse(x, fmm->tree);
	return val;
}


double RBF::fmmComputeValueRecurse(vec3 x, BHNode *myNode)
{
	double val=0;
  
	if(length(myNode->center - x)/length(myNode->box.min-myNode->box.max) < 0)  //far away
	{
		val = myNode->coeff*fmmComputeKernel(x, myNode);
	}
	else if (!myNode->leaf) // leaf
	{
		for(int i=0; i<8; i++)
		{
			if(myNode->nodes[i] != 0)
			{
				val+=fmmComputeValueRecurse(x, myNode->nodes[i]);
			}
		}
	}
	else //close, but not a leaf
	{
		int n=myNode->pts.size();
		for(int i=0; i<n; i++)
		{
			val+=coeff[myNode->pts[i]]*computeKernel(myNode->pts[i], x);
		}
	}
	return val;
}

double RBF::fmmComputeKernel(vec3 b, BHNode *myNode)
{
	double r = length(myNode->center - b);
  
	return computeRadialFunction(r);
}


double RBF::fmmComputeKernel(BHNode *myNode, vec3 b)
{
	return fmmComputeKernel(b, myNode);
}
