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
#include <rbf/Application/Tools/src/Surface.h>
#include <rbf/Application/Tools/src/RBF.h>
#include <rbf/Application/Tools/src/ScatteredData.h>

#include <rbf/Core/Isosurface/Isosurface.h>

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
using namespace Plugin::Core;

// TODO really worth having a private class here?
class ActionRBFPrivate
{
public:
  std::string targetLayerID_;
  VertexList vertices_;
  double isovalue_;
  bool capIsosurface_;
  std::string kernel_;
  double x_, y_, z_;
};

class RBFAlgo : public LayerFilter
{
  
public:
  LayerHandle srcLayer_;
  LayerHandle dstLayer_;
  ActionRBFPrivateHandle actionInternal_;
  Point origin_, min_, max_;
  double dataRange_;

  void augmentNormalData(ScatteredData *data, double vectorLength);
  vec3 findNormal(ScatteredData *data, int n);
  void computeRange();

  static const int DIMS;

	RBFAlgo();
	virtual ~RBFAlgo();
  
  SCI_BEGIN_RUN()
  {
    this->computeRange();
    DataLayerHandle srcDataLayer = boost::dynamic_pointer_cast<DataLayer>(srcLayer_);
    DataLayerHandle dstDataLayer = boost::dynamic_pointer_cast<DataLayer>(dstLayer_);
    GridTransform gridTransform = srcDataLayer->get_grid_transform();
    DataBlockHandle data = srcDataLayer->get_data_volume()->get_data_block();

//    if (! this->create_and_lock_data_layer_test( dst_layer, "surface", gridTransform, meta_data ) )
//    {
//      this->report_error( "Could not create layer." );
//    }
    DataVolumeHandle volume;
    if ( !( DataVolume::DuplicateVolume( srcDataLayer->get_data_volume(), volume ) ) )
    {
      this->report_error( "Could not duplicate volume." );
      return;
    }
//    DataVolumeHandle volume = DataVolumeHandle( new DataVolume( gridTransform, dstDataBlock ) );
    if (! this->dispatch_insert_data_volume_into_layer( this->dstLayer_, volume, true ) )
    {
      this->report_error( "Could not insert layer." );
      return;
    }

    double spacing_x = (this->max_.x() - this->min_.x()) / static_cast<double>(actionInternal_->x_),
           spacing_y = (this->max_.y() - this->min_.y()) / static_cast<double>(actionInternal_->y_),
           spacing_z = (this->max_.z() - this->min_.z()) / static_cast<double>(actionInternal_->z_);
    std::cerr << "spacing=" << spacing_x << "," << spacing_y << "," << spacing_z << std::endl;
    std::cerr << "origin=" << this->origin_ << std::endl;

    // temporary: destination layer will be created in the usual way
//    LayerHandle dst_layer;
//    LayerMetaData meta_data;
//    Vector test_spacing(0.6, 0.5, 0.1);
//    Vector dims(100, 100, 100);
//    Point origin(-30, -50, 80);
//    Transform transform(origin,
//                        Vector( test_spacing.x(), 0.0 , 0.0 ),
//                        Vector( 0.0, test_spacing.y(), 0.0 ),
//                        Vector( 0.0, 0.0, test_spacing.z() ));
//    GridTransform surfaceGridTransform( dims.x(), dims.y(), dims.z(), transform );
//    surfaceGridTransform.set_originally_node_centered( false );

    //Take the input
    ScatteredData *mySurfaceData = new ScatteredData();

    //	readSurfaceDataFile("./segmentation/sample_points_DEMRI.txt", mySurfaceData);
    std::vector<Surface*> mySurface;

    for (int j = 0; j < this->actionInternal_->vertices_.size(); j++)
    {
//      for(int i = 0; i < DIMS; i++)
//        mySurfaceData->x[i].push_back(myData[j][i]);

      mySurfaceData->x[0].push_back( this->actionInternal_->vertices_[j].x() );
      mySurfaceData->x[1].push_back( this->actionInternal_->vertices_[j].y() );
      mySurfaceData->x[2].push_back( this->actionInternal_->vertices_[j].z() );

      // TODO: hardcoded isovalue
      mySurfaceData->fnc.push_back(0);
    }

    // TODO: ???
//    for(int i = 0; i < DIMS; i++)
//      mySurfaceData->x[i].pop_back();

    mySurfaceData->fnc.pop_back();

//    double dataRange = 20;
    // vector length ~ 5% of data range
    //double vectorLength = dataRange*0.05;
//    std::cerr << "data range=" << dataRange << ", vector length=" << vectorLength << std::endl;
//    double vectorLength = 10;
    double vectorLength = 1;
    augmentNormalData(mySurfaceData, vectorLength);
    
    RBF *mySurfaceRBF;
    Kernel myKernel = ThinPlate;
    
    mySurfaceRBF = new RBF(mySurfaceData, myKernel);
    mySurfaceRBF->setDataReduction(Random);
    
    mySurface.push_back(new Surface(mySurfaceData, mySurfaceRBF));
    
    
    //Construct RBFs
    mySurface[0]->computeRBF();

    //Compute material
    std::vector<std::vector<std::vector<double> > > value;
    value.resize(this->actionInternal_->x_);
//    vec3 myOrigin(origin.x(), origin.y(), origin.z());
    vec3 myOrigin(this->origin_.x(), this->origin_.y(), this->origin_.z());
    for (int i = 0; i < this->actionInternal_->x_; i++)
    {
//      if (i % 10 == 0)
//        printf("%d/100 done\n", i); fflush(stdout);
      value[i].resize(this->actionInternal_->y_);
      for(int j = 0; j < this->actionInternal_->y_; j++)
      {
        //if(j%10==0)
        //	printf("\t%d/100 done\n", j); fflush(stdout);
        value[i][j].resize(this->actionInternal_->z_);
        for(int k = 0; k < this->actionInternal_->z_; k++)
        {
          //if(k%10==0)
          //	printf("\t\t%d/100 done\n", k); fflush(stdout);
          vec3 location = myOrigin + spacing_x * i * vec3::unitX + spacing_y * j * vec3::unitY + spacing_z * k * vec3::unitZ;
//          vec3 location = myOrigin + test_spacing.x() * i * vec3::unitX + test_spacing.y() * j * vec3::unitY + test_spacing.z() * k * vec3::unitZ;
          //std::cout<<"Computing Val ... "<<std::endl;
          double myVal = mySurface[0]->computeValue(location);
          //printf("Interpolant: %lf %lf %lf %lf\n", location[0], location[1], location[2], myVal); fflush(stdout);
          value[i][j][k]=myVal;
        }
      }
    }

    GridTransform surfaceGridTransform( actionInternal_->x_, actionInternal_->y_, actionInternal_->z_, gridTransform.transform() );
    surfaceGridTransform.set_originally_node_centered( gridTransform.get_originally_node_centered() );
    DataBlockHandle surfaceDataBlock = StdDataBlock::New( this->actionInternal_->x_, this->actionInternal_->y_, this->actionInternal_->z_, DataType::FLOAT_E );
    DataVolumeHandle surfaceDataVolume = DataVolumeHandle( new DataVolume( surfaceGridTransform, surfaceDataBlock ) );

    // Initialize the result data with 0
    float* surfaceData = reinterpret_cast< float* >( surfaceDataBlock->get_data() );
    size_t total_voxels = surfaceDataBlock->get_size();
    memset( surfaceData, 0, total_voxels );

    for (size_t x = 0; x < this->actionInternal_->x_; x++)
    {
      for (size_t y = 0; y < this->actionInternal_->y_; y++)
      {
        for (size_t z = 0; z < this->actionInternal_->z_; z++)
        {
          size_t index = surfaceDataBlock->to_index( x, y, z );
          surfaceData[index] = value[x][y][z];
          //printf("%d %d %d %lf\n", x,y,z,value[x][y][z]);
        }
      }
    }

//    if (! this->dispatch_unlock_layer( this->dstLayer_ ) )
//    {
//      this->report_error( "Could not unlock layer." );
//    }

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
      nrrd_file << "sizes: " << this->actionInternal_->x_ << " " << this->actionInternal_->y_ << " " << this->actionInternal_->z_ << std::endl;
      nrrd_file << "axis mins: " << 0 << ", " << 0 << ", " << 0 << std::endl;
//      nrrd_file << "spacings: " << test_spacing.x() << " " << test_spacing.y() << " " << test_spacing.z() << std::endl;
      nrrd_file << "spacings: " << spacing_x << " " << spacing_y << " " << spacing_z << std::endl;
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
            float val = value[i][j][k];
            nrrd_file.write((char*)&val, sizeof(float));
          }
        }
      }
      nrrd_file.close();
    }

    IsosurfaceHandle iso = IsosurfaceHandle( new Isosurface(surfaceDataVolume, this->actionInternal_->isovalue_) );
    iso->compute( 1.0, this->actionInternal_->capIsosurface_, boost::bind( &LayerFilter::check_abort, this ) );
    dstDataLayer->set_isosurface(iso);
    iso->export_stl_isosurface("surface.stl", "test_surface");
    iso->export_vtk_isosurface("surface.vtk");
    std::cerr << "Isosurface surface area: " << iso->surface_area() << std::endl;
    if (iso->surface_area() == 0)
    {
      this->report_error( "Empty isosurface." );
      return;
    }

    ViewerHandle viewer = ViewerManager::Instance()->get_viewer( "viewer0" );
    if (! viewer)
    {
      this->report_error( "Bad viewer handle." );
      return;
    }
    ActionSet::Dispatch( Interface::GetWidgetActionContext(), viewer->volume_isosurfaces_visible_state_, true );

    ActionSet::Dispatch( Interface::GetWidgetActionContext(), dstDataLayer->show_isosurface_state_, true );
    ActionSet::Dispatch( Interface::GetWidgetActionContext(), dstDataLayer->iso_generated_state_, true );
    dstDataLayer->isosurface_updated_signal_();
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

RBFAlgo::RBFAlgo() : dataRange_(-1)
{
}

RBFAlgo::~RBFAlgo()
{
}

const int RBFAlgo::DIMS = 3;

void RBFAlgo::computeRange()
{
  VertexList::iterator it = this->actionInternal_->vertices_.begin();
  double minX = it->x(), minY = it->y(), minZ = it->z();
  double maxX = it->x(), maxY = it->y(), maxZ = it->z();
  ++it;

  for (; it != this->actionInternal_->vertices_.end(); ++it)
  {
    minX = Min(minX, it->x());
    maxX = Max(maxX, it->x());
    minY = Min(minY, it->y());
    maxY = Max(maxY, it->y());
    minZ = Min(minZ, it->z());
    maxZ = Max(maxZ, it->z());
  }

  std::cerr << "min: " << minX << " " << minY << " " << minZ << std::endl;
  std::cerr << "max: " << maxX << " " << maxY << " " << maxZ << std::endl;

  this->min_.x(minX);
  this->min_.y(minY);
  this->min_.z(minZ);
  this->max_.x(maxX);
  this->max_.y(maxY);
  this->max_.z(maxZ);

  double dataRangeX = maxX - minX;
  double dataRangeY = maxY - minY;
  double dataRangeZ = maxZ - minZ;
  this->dataRange_ = (dataRangeX + dataRangeY + dataRangeZ)/3;

  std::cerr << "dataRange: " << dataRangeX << " " << dataRangeY << " " << dataRangeZ << ", " << this->dataRange_ << std::endl;

  this->origin_.x((minX + maxX)/2);
  this->origin_.y((minY + maxY)/2);
  this->origin_.z((minZ + maxZ)/2);
  std::cerr << "origin: " << this->origin_ << std::endl;
}

vec3 RBFAlgo::findNormal(ScatteredData *data, int n)
{
  int tot = data->x[0].size();
  int prev = (n-1) >= 0 ? n-1 : tot-1;
  int next = (n+1) < tot ? n+1 : 0;

  while( data->x[2][prev] != data->x[2][n] )
  {
    prev = (prev-1) >= 0 ? prev-1 : tot-1;
  }

  while( data->x[2][next] != data->x[2][n] )
  {
    next = (next+1) < tot ? next+1 : 0;
  }
//  printf("%d %d %d %d\n", prev,n,next,tot); fflush(stdout);

  vec3 a(data->x[0][n], data->x[1][n], data->x[2][n]);
  vec3 b(data->x[0][prev], data->x[1][prev], data->x[2][prev]);
  vec3 c(data->x[0][next], data->x[1][next], data->x[2][next]);
  vec3 one = b-a;
  vec3 two = c-a;
  vec3 ret = one+two;

  return ret;
}

void RBFAlgo::augmentNormalData(ScatteredData *data, double vectorLength)
{
  const int N = data->x[0].size();
  for(int i = 0; i < N; i++)
  {
    vec3 myNormal = findNormal(data, i);
    myNormal = normalize(myNormal);

    for(int j = 0; j < DIMS; j++)
    {
      data->x[j].push_back(data->x[j][i] + myNormal[j]);
    }

//    data->fnc.push_back(10);
    data->fnc.push_back(vectorLength);

    for(int j = 0; j < DIMS; j++)
    {
      data->x[j].push_back(data->x[j][i] - myNormal[j]);
    }
    // TODO: value of function one unit away - needs to be generalized!
//    data->fnc.push_back(-10);
    data->fnc.push_back(-vectorLength);
  }
}

ActionRBF::ActionRBF() :
  private_( new ActionRBFPrivate )
{
  this->add_layer_id( this->private_->targetLayerID_ );
  this->add_parameter( this->private_->vertices_ );
  this->add_parameter( this->private_->isovalue_ );
  this->add_parameter( this->private_->kernel_ );
  this->add_parameter( this->private_->capIsosurface_ );
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
                         double isovalue,
                         const std::string& kernel,
                         bool cap_isosurface,
                         double x,
                         double y,
                         double z
                        )
{
  ActionRBF* action = new ActionRBF;
  action->private_->targetLayerID_ = target;
  action->private_->vertices_ = vertices;
  action->private_->isovalue_ = isovalue;
  action->private_->kernel_ = kernel;
  action->private_->capIsosurface_ = cap_isosurface;
  action->private_->x_ = x;
  action->private_->y_ = y;
  action->private_->z_ = z;

  ActionDispatcher::PostAction( ActionHandle( action ), context );
}

}}
