/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2014 Scientific Computing and Imaging Institute,
 University of Utah.
 
 
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
*/

#include <Core/Action/ActionFactory.h>

//#include <Application/Provenance/Provenance.h>
//#include <Application/ToolManager/ToolManager.h>
#include <Application/Filters/LayerFilter.h>
#include <Core/Utils/Log.h>
#include <Core/Math/MathFunctions.h>

#include <rbf/Application/Tools/Actions/ActionRBF.h>

// RBF library includes
#include <rbf/Application/Tools/src/RBFInterface.h>
#include <rbf/Application/Tools/src/vec3.h>

//#include <Core/Isosurface/Isosurface.h>

// test
#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <Core/DataBlock/StdDataBlock.h>
#include <Application/ViewerManager/ViewerManager.h>
#include <Core/State/Actions/ActionSetAt.h>
#include <Core/Interface/Interface.h>
#include <Core/Volume/DataVolume.h>
// test

CORE_REGISTER_ACTION( Plugin::Application, RBF )

namespace Plugin
{

namespace Application
{

using namespace ::Seg3D;
using namespace ::Core;

// TODO really worth having a private class here?
class ActionRBFPrivate
{
public:
  std::string targetLayerID_;
  VertexList vertices_;
//  double isovalue_;
//  bool capIsosurface_;
  std::string kernel_;
  double x_, y_, z_;
};

class RBFAlgo : public LayerFilter
{
  
public:
  LayerHandle srcLayer_;
  LayerHandle dstLayer_;
  ActionRBFPrivateHandle actionInternal_;
//  Point origin_, min_, max_;
//  double dataRange_;
//
//  void augmentNormalData(ScatteredData *data, double vectorLength);
//  vec3 findNormal(ScatteredData *data, int n);
//  void computeRange();
//
//  static const int DIMS;

	RBFAlgo();
	virtual ~RBFAlgo();
  
  SCI_BEGIN_RUN()
  {
    DataLayerHandle srcDataLayer = boost::dynamic_pointer_cast<DataLayer>(srcLayer_);
    DataLayerHandle dstDataLayer = boost::dynamic_pointer_cast<DataLayer>(dstLayer_);
    GridTransform gridTransform = srcDataLayer->get_grid_transform();

    std::vector<vec3> rbfPointData;
    for (VertexList::iterator it = this->actionInternal_->vertices_.begin();
         it != this->actionInternal_->vertices_.end();
         ++it)
    {
      rbfPointData.push_back( vec3(it->x(), it->y(), it->z() ) );
    }
    // origin and size from source data layer
    Point origin = gridTransform.origin();
    vec3 rbfOrigin(origin.x(), origin.y(), origin.z());
    vec3 rbfGridSize(gridTransform.get_nx(), gridTransform.get_ny(), gridTransform.get_nz());
    vec3 rbfSampling(this->actionInternal_->x_, this->actionInternal_->y_, this->actionInternal_->z_);

    RBFInterface rbfAlgo(rbfPointData, rbfOrigin, rbfGridSize, rbfSampling);
    DataStructure rbfData = rbfAlgo.value;
    Core::DataBlockHandle dstDataBlock = Core::StdDataBlock::New( gridTransform, Core::DataType::FLOAT_E );
    if ( ! dstDataBlock )
    {
      this->report_error( "Could not allocate enough memory." );
      return;
    }

//    // test
//    for (size_t i = 0; i < dstDataBlock->get_nx(); ++i)
//    {
//      for (size_t j = 0; j < dstDataBlock->get_ny(); ++j)
//      {
//        for (size_t k = 0; k < dstDataBlock->get_nz(); ++k)
//        {
//          dstDataBlock->set_data_at( i, j, k, 1 );
//        }
//      }
//    }
//    // test

    for (size_t i = 0; i < dstDataBlock->get_nx(); ++i)
    {
      for (size_t j = 0; j < dstDataBlock->get_ny(); ++j)
      {
        for (size_t k = 0; k < dstDataBlock->get_nz(); ++k)
        {
          // dims do not appear to match - this should not be necessary!
          if (i >= rbfAlgo.nx || j >= rbfAlgo.ny || k >= rbfAlgo.nz)
          {
            // test
//            dstDataBlock->set_data_at( i, j, k, 1 );
            continue;
          }
          dstDataBlock->set_data_at( i, j, k, rbfData[i][j][k] );
        }
      }
    }

    this->dispatch_insert_data_volume_into_layer(
      this->dstLayer_,
      Core::DataVolumeHandle(new Core::DataVolume( this->dstLayer_->get_grid_transform(), dstDataBlock ) ),
      true );

    // test
    std::string filename = "surface.nrrd";
    std::cout << "Writing file '" << filename << "'" << std::endl;
    std::ofstream nrrd_file(filename.c_str(), std::ofstream::binary);
    
    if (nrrd_file.is_open())
    {
      nrrd_file << "NRRD0001" << std::endl;
      nrrd_file << "# Complete NRRD file format specification at:" << std::endl;
      nrrd_file << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
      nrrd_file << "type: float" << std::endl;
      nrrd_file << "dimension: 3" << std::endl;
      nrrd_file << "sizes: " << rbfAlgo.nx << " " << rbfAlgo.ny << " " << rbfAlgo.nz << std::endl;
      nrrd_file << "axis mins: " << 0 << ", " << 0 << ", " << 0 << std::endl;
      nrrd_file << "spacings: " << rbfAlgo.spacing_x << " " << rbfAlgo.spacing_y << " " << rbfAlgo.spacing_z << std::endl;
      nrrd_file << "centerings: cell cell cell" << std::endl;
      nrrd_file << "endian: little" << std::endl;
      nrrd_file << "encoding: raw" << std::endl;
      nrrd_file << std::endl;
      
      // write data portion
      for(int k=0; k < this->actionInternal_->z_; k++)
      {
        for(int j=0; j < this->actionInternal_->y_; j++)
        {
          for(int i=0; i < this->actionInternal_->x_; i++)
          {
            float val = rbfAlgo.value[i][j][k];
            nrrd_file.write((char*)&val, sizeof(float));
          }
        }
      }
      nrrd_file.close();
    }
    // test
}
  SCI_END_RUN()
  
  // GET_FITLER_NAME:
  // The name of the filter, this information is used for generating new layer labels.
  virtual std::string get_filter_name() const
  {
    return "RBF Tool";
  }
  
  // GET_LAYER_PREFIX:
  // This function returns the name of the filter. The latter is prepended to the new layer name,
  // when a new layer is generated.
  virtual std::string get_layer_prefix() const
  {
    return "RBF";	
  }	
};

RBFAlgo::RBFAlgo()// : dataRange_(-1)
{
}

RBFAlgo::~RBFAlgo()
{
}

//const int RBFAlgo::DIMS = 3;

//void RBFAlgo::computeRange()
//{
//  VertexList::iterator it = this->actionInternal_->vertices_.begin();
//  double minX = it->x(), minY = it->y(), minZ = it->z();
//  double maxX = it->x(), maxY = it->y(), maxZ = it->z();
//  ++it;
//
//  for (; it != this->actionInternal_->vertices_.end(); ++it)
//  {
//    minX = Min(minX, it->x());
//    maxX = Max(maxX, it->x());
//    minY = Min(minY, it->y());
//    maxY = Max(maxY, it->y());
//    minZ = Min(minZ, it->z());
//    maxZ = Max(maxZ, it->z());
//  }
//
//  std::cerr << "min: " << minX << " " << minY << " " << minZ << std::endl;
//  std::cerr << "max: " << maxX << " " << maxY << " " << maxZ << std::endl;
//
//  this->min_.x(minX);
//  this->min_.y(minY);
//  this->min_.z(minZ);
//  this->max_.x(maxX);
//  this->max_.y(maxY);
//  this->max_.z(maxZ);
//
//  double dataRangeX = maxX - minX;
//  double dataRangeY = maxY - minY;
//  double dataRangeZ = maxZ - minZ;
//  this->dataRange_ = (dataRangeX + dataRangeY + dataRangeZ)/3;
//
//  std::cerr << "dataRange: " << dataRangeX << " " << dataRangeY << " " << dataRangeZ << ", " << this->dataRange_ << std::endl;
//
//  this->origin_.x((minX + maxX)/2);
//  this->origin_.y((minY + maxY)/2);
//  this->origin_.z((minZ + maxZ)/2);
//  std::cerr << "origin: " << this->origin_ << std::endl;
//}
//
//vec3 RBFAlgo::findNormal(ScatteredData *data, int n)
//{
//  int tot = data->x[0].size();
//  int prev = (n-1) >= 0 ? n-1 : tot-1;
//  int next = (n+1) < tot ? n+1 : 0;
//
//  while( data->x[2][prev] != data->x[2][n] )
//  {
//    prev = (prev-1) >= 0 ? prev-1 : tot-1;
//  }
//
//  while( data->x[2][next] != data->x[2][n] )
//  {
//    next = (next+1) < tot ? next+1 : 0;
//  }
////  printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);
//
//  vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
//  vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
//  vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
//  vec3 one = b-a;
//  vec3 two = c-a;
//  vec3 ret = one+two;
//
//  return ret;
//}
//
//void RBFAlgo::augmentNormalData(ScatteredData *data, double vectorLength)
//{
//  const int N = data->x[0].size();
//  for(int i = 0; i < N; i++)
//  {
//    vec3 myNormal = findNormal(data, i);
//    myNormal = normalize(myNormal);
//
//    for(int j = 0; j < DIMS; j++)
//    {
//      data->x[j].push_back(data->x[j][i] + myNormal[j]);
//    }
//
////    data->fnc.push_back(10);
//    data->fnc.push_back(vectorLength);
//
//    for(int j = 0; j < DIMS; j++)
//    {
//      data->x[j].push_back(data->x[j][i] - myNormal[j]);
//    }
//    // TODO: value of function one unit away - needs to be generalized!
////    data->fnc.push_back(-10);
//    data->fnc.push_back(-vectorLength);
//  }
//}

ActionRBF::ActionRBF() :
  private_( new ActionRBFPrivate )
{
  this->add_layer_id( this->private_->targetLayerID_ );
  this->add_parameter( this->private_->vertices_ );
//  this->add_parameter( this->private_->isovalue_ );
  this->add_parameter( this->private_->kernel_ );
//  this->add_parameter( this->private_->capIsosurface_ );
  this->add_parameter( this->private_->x_ );
  this->add_parameter( this->private_->y_ );
  this->add_parameter( this->private_->z_ );
  this->add_parameter( this->sandbox_ );
}

bool ActionRBF::validate( ActionContextHandle& context )
{
  if (this->private_->vertices_.size() < 2)
  {
    context->report_error("Non-trivial number of points needed.");
    return false;
  }
  return true;
}

bool ActionRBF::run( ActionContextHandle& context, ActionResultHandle& result )
{
	boost::shared_ptr< RBFAlgo > algo( new RBFAlgo() );
  
	// Set up parameters
	algo->set_sandbox( this->sandbox_ );
  algo->actionInternal_ = private_;

	// Find the handle to the layer
	if ( !( algo->find_layer( this->private_->targetLayerID_, algo->srcLayer_ ) ) )
	{
		return false;
	}

  // no replace option

  // Lock the src layer, so it cannot be used else where
  algo->lock_for_use( algo->srcLayer_ );

  // Create the destination layer, which will show progress
  algo->create_and_lock_data_layer_from_layer( algo->srcLayer_, algo->dstLayer_ );

	// Return the id of the destination layer.
	result = ActionResultHandle( new ActionResult( algo->dstLayer_->get_layer_id() ) );

	// If the action is run from a script (provenance is a special case of script),
	// return a notifier that the script engine can wait on.
	if ( context->source() == ActionSource::SCRIPT_E ||
       context->source() == ActionSource::PROVENANCE_E )
	{
		context->report_need_resource( algo->get_notifier() );
	}

	// Build the undo-redo record
	algo->create_undo_redo_and_provenance_record( context, this->shared_from_this() );

	// Start the filter.
	Runnable::Start( algo );
  
	return true;
}

void ActionRBF::Dispatch(
                         ActionContextHandle context,
                         const std::string& target,
                         const VertexList& vertices,
//                         double isovalue,
                         const std::string& kernel,
//                         bool cap_isosurface,
                         double x,
                         double y,
                         double z
                        )
{
  ActionRBF* action = new ActionRBF;
  action->private_->targetLayerID_ = target;
  action->private_->vertices_ = vertices;
//  action->private_->isovalue_ = isovalue;
  action->private_->kernel_ = kernel;
//  action->private_->capIsosurface_ = cap_isosurface;
  action->private_->x_ = x;
  action->private_->y_ = y;
  action->private_->z_ = z;

  ActionDispatcher::PostAction( ActionHandle( action ), context );
}

}}
