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

#ifndef _RBF_H_
#define _RBF_H_

#include <rbf/Application/Tools/src/ScatteredData.h>
#include <rbf/Application/Tools/src/SparseMatrix.h>
#include <rbf/Application/Tools/src/vec3.h>
#include <rbf/Application/Tools/src/FMM.h>

// STL Includes
#include <vector>
#include <utility>

enum Kernel {Guassian, ThinPlate, MultiQuadratic};
enum Acceleration {None, FastMultipole};
enum DataReduction {All, Random};

class RBF
{
public:
	RBF();
	RBF(ScatteredData *myData, Kernel myKernel);
  
	void setKernel(Kernel myKernel);
	void setData(ScatteredData *myData);
	void setAcceleration(Acceleration myAcceleration);
	void setDataReduction(DataReduction myDataReduction);
  
	void computeFunction();
	double computeValue(vec3 x);

private:
	ScatteredData *data, *completeData;
	Kernel kernel;
  std::vector<double> coeff;
	Acceleration acceleration;
	DataReduction dataReduction;
  
	FMM *fmm;
  
	void computeFunctionForData();
	void computeErrorForData(std::vector<std::pair<double, int> > &error);
  
	double computeKernel(int i, int j);
	double computeKernel(int i, vec3 b);
	double computeRadialFunction(double r);
  
  
	void fmmComputeFunction();
	void fmmBuildTree();
	void fmmPrintTree(BHNode* myNode, int stack);
	void fmmBuildTree(std::vector<int> &myPoints, BHNode *myNode);
	double fmmComputeValue(vec3 x);
	double fmmComputeValueRecurse(vec3 x, BHNode *myNode);
	double fmmComputeKernel(vec3 b, BHNode *myNode);
	double fmmComputeKernel(BHNode *myNode, vec3 b);
	
};

#endif //_RBF_H_
