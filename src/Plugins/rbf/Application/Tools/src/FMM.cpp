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

FMM::FMM()
{
	theta=0.75;
}

FMM::~FMM()
{
	freeTheTree(tree);
}

void FMM::freeTheTree(BHNode *myNode)
{
	if (myNode == 0)
		return;
  
	for (int i=0; i<8; i++)
		freeTheTree(myNode->nodes[i]);
}

void FMM::setTheta(double myTheta)
{
	theta = myTheta;
}
