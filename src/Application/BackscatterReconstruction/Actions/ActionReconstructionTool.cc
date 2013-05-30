/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#include <Core/Application/Application.h>
#include <Core/Action/ActionDispatcher.h>
#include <Core/Action/ActionFactory.h>
#include <Core/State/Actions/ActionSet.h>

// Application Includes
#include <Application/BackscatterReconstruction/Actions/ActionReconstructionTool.h>
#include <Application/BackscatterReconstruction/Algorithm/recon_api.h>

#include <Application/ToolManager/Actions/ActionOpenTool.h>
#include <Application/ViewerManager/ViewerManager.h>
#include <Application/Filters/ITKFilter.h>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/conversion.hpp>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, ReconstructionTool )

// TODO: need more systematic way to check number of layers
#define LAYER_COUNT 2

namespace Seg3D
{

class ReconstructionProgressRunnable : public Core::Runnable
{
  ActionReconstructionTool::callback_type callback_;
  long msTimeInterval_;
  bool stop_;

public:
  ReconstructionProgressRunnable(ActionReconstructionTool::callback_type callback, long msTimeInterval=500)
    : callback_(callback), msTimeInterval_(msTimeInterval), stop_(false) {}
  
  void stop() { stop_ = true; }

  virtual void run()
  {
    // TODO: temporary? See TODO below...
    double scale = 1000;
    while (! stop_ )
    {
      ITKFilter::UCHAR_IMAGE_TYPE::Pointer reconVolum = ITKFilter::UCHAR_IMAGE_TYPE::New();
      double progress = ReconstructionGetMaterialVolume(reconVolum);
      if ( this->callback_ != 0 )
      {
        // TODO: using scale temporarily because progress values are < 0.01
        // before reconstruction either exits or crashes
        Core::Application::PostEvent( boost::bind( this->callback_, progress*scale, 0, 1.0 ) );
      }
      // test
//      ITKFilter::UCHAR_IMAGE_TYPE::SizeType sizes = reconVolum->GetLargestPossibleRegion().GetSize();
//      const size_t SIZE = sizes[0] * sizes[1] * sizes[2];
//      for (size_t i = 0; i < SIZE; ++i)
//      {
//        if (reconVolum->GetBufferPointer()[i] != 0)
//        {
//          printf("reconVolum[%lu]=%hhu\n", i, reconVolum->GetBufferPointer()[i]);
//        }
//      }

      boost::this_thread::sleep(boost::posix_time::milliseconds(msTimeInterval_));
    }
    // make sure progress bar completes
    Core::Application::PostEvent( boost::bind( this->callback_, 100.0, 0, 1.0 ) );
  }
};


class ReconstructionToolAlgo : public ITKFilter
{  

public:
  LayerHandle src_layer_;
  std::vector< LayerHandle > initialGuessSet_;

  std::vector<LayerHandle> dst_layer_;

  std::string outputDir_;
  double xyVoxelSizeScale_;
  double zVoxelSizeScale_;
  int iterations_;
  ActionReconstructionTool::callback_type callback_;
  
public:
  // RUN:
  // Implementation of run of the Runnable base class, this function is called when the thread
  // is launched.
  
  // NOTE: The macro needs a data type to select which version to run. This needs to be
  // a member variable of the algorithm class.
  SCI_BEGIN_ITK_RUN()
  {
    boost::filesystem::path algorithm_work_dir;
    boost::filesystem::path algorithm_config_file;    
    boost::filesystem::path algorithm_geometry_file;
    Core::Application::Instance()->get_algorithm_config(algorithm_work_dir,
                                                        algorithm_config_file,
                                                        algorithm_geometry_file);
    std::cerr << "work dir=" << algorithm_work_dir.string() << ", config file="
      << algorithm_config_file << ", geometry_file=" << algorithm_geometry_file << std::endl;
    
    Core::ITKImageDataT<float>::Handle image; 
    this->get_itk_image_from_layer<float>( this->src_layer_, image );
    ITKFilter::FLOAT_IMAGE_TYPE::SizeType inSize = image->get_image()->GetLargestPossibleRegion().GetSize();

    DataLayerHandle data_layer = boost::dynamic_pointer_cast<DataLayer>( this->src_layer_ );
    Core::VolumeHandle srcVolumeHandle = data_layer->get_volume();
    
    ITKFilter::UCHAR_IMAGE_TYPE::Pointer initialGuess;

    if (this->initialGuessSet_.size() > 0)
    {
      initialGuess = ITKFilter::UCHAR_IMAGE_TYPE::New();
      ITKFilter::UCHAR_IMAGE_TYPE::SizeType size;
      size.SetElement(0, inSize[0]);
      size.SetElement(1, inSize[1]);
      size.SetElement(2, inSize[2]);
      initialGuess->SetRegions(ITKFilter::UCHAR_IMAGE_TYPE::RegionType(size));
      initialGuess->Allocate();

      // TODO: ordering of mask labels?
      for (size_t i = 0; i < this->initialGuessSet_.size(); ++i)
      {
        Core::MaskLayerHandle mask_layer = boost::dynamic_pointer_cast<MaskLayer>( this->initialGuessSet_[i] );
        Core::MaskVolumeHandle volumeHandle = mask_layer->get_mask_volume();
        if ( volumeHandle->get_size() != srcVolumeHandle->get_size() )
        {
          this->report_error("Mask size does not match data volume size");
          return;
        }
        Core::MaskDataBlockHandle dataBlock = volumeHandle->get_mask_data_block();
        
        for (size_t j = 0; j < volumeHandle->get_size(); j++)
        {
          if (dataBlock->get_mask_at(j))
          {
            // set disk IDs from 1 onwards
            initialGuess->GetBufferPointer()[j] = (i+1);
          }
          else
          {
            initialGuess->GetBufferPointer()[j] = 0;
          }
        }
      }
    }

    // allocated by algorithm
    ITKFilter::UCHAR_IMAGE_TYPE::Pointer finalReconVolume;

    float voxelSizeCM[3] = { this->xyVoxelSizeScale_, this->xyVoxelSizeScale_, this->zVoxelSizeScale_ };
    
    boost::shared_ptr<ReconstructionProgressRunnable> progress(
      new ReconstructionProgressRunnable(this->callback_) );
    Core::Runnable::Start( progress );

    ReconstructionStart(image->get_image(),
                        initialGuess,
                        algorithm_config_file.string().c_str(),
                        voxelSizeCM,
                        iterations_,
                        finalReconVolume);
    progress->stop();

    
    for (size_t i = 0; i < LAYER_COUNT; ++i)
    {
      MaskLayerHandle mask_layer = boost::dynamic_pointer_cast<MaskLayer>( this->dst_layer_[ i ] );
      Core::MaskDataBlockHandle maskDataBlock;
      Core::MaskDataBlockManager::Create( this->dst_layer_[ i ]->get_grid_transform(), maskDataBlock );
      Core::MaskVolumeHandle mask_volume( new Core::MaskVolume(this->dst_layer_[ i ]->get_grid_transform(), maskDataBlock ) );
      
      for (size_t j = 0; j < maskDataBlock->get_size(); j++ )
      {
        if ( finalReconVolume->GetBufferPointer()[j] == (i+1) )
        {
          maskDataBlock->set_mask_at(j);
        }
      }
      
      if (! this->dispatch_insert_mask_volume_into_layer( this->dst_layer_[ i ], mask_volume ) )
      {
        std::ostringstream oss;
        oss << "Could not insert mask volume " << this->dst_layer_[ i ]->get_layer_id() << " into layer";
        this->report_error(oss.str());
      }
    }

  }
  SCI_END_ITK_RUN()
  
  // GET_FITLER_NAME:
  // The name of the filter, this information is used for generating new layer labels.
  virtual std::string get_filter_name() const
  {
    return "ReconstructionToolAlgorithm";
  }
  
  // GET_LAYER_PREFIX:
  // This function returns the name of the filter. The latter is prepended to the new layer name, 
  // when a new layer is generated. 
  virtual std::string get_layer_prefix() const
  {
    return "ReconstructionTool";
  }
};

bool ActionReconstructionTool::validate( Core::ActionContextHandle& context )
{
  // Make sure that the sandbox exists
  if ( !LayerManager::CheckSandboxExistence( this->sandbox_, context ) )
  {
    return false;
  }
  
  // Check for layer existence and type information
  if ( ! LayerManager::CheckLayerExistenceAndType( this->target_layer_, 
    Core::VolumeType::DATA_E, context, this->sandbox_ ) ) return false;
  
  // Check for layer availability 
  if ( ! LayerManager::CheckLayerAvailabilityForProcessing( this->target_layer_, 
    context, this->sandbox_ ) ) return false;
  
  return true; // validated
}

ActionReconstructionTool::ActionReconstructionTool()
{
  this->add_layer_id( this->target_layer_ );
  this->add_parameter( this->initalGuessSet_ );
  this->add_parameter( this->iterations_ );
  this->add_parameter( this->outputDir_ );
  this->add_parameter( this->xyVoxelSizeScale_ );
  this->add_parameter( this->zVoxelSizeScale_ );
  this->add_parameter( this->sandbox_ );
}

bool ActionReconstructionTool::run( Core::ActionContextHandle& context,
                                    Core::ActionResultHandle& result )
{
  // Create algorithm
  boost::shared_ptr<ReconstructionToolAlgo> algo( new ReconstructionToolAlgo );

  // Copy the parameters over to the algorithm that runs the filter
  algo->set_sandbox( this->sandbox_ );

  // Set up parameters
  algo->outputDir_ = this->outputDir_;
  algo->xyVoxelSizeScale_ = this->xyVoxelSizeScale_;
  algo->zVoxelSizeScale_ = this->zVoxelSizeScale_;
  algo->iterations_ = this->iterations_;
  algo->callback_ = this->callback_;

  // Find the handle to the layer
  algo->src_layer_ = LayerManager::FindLayer( this->target_layer_, this->sandbox_ );
  // Check if layer really exists
  if ( ! algo->src_layer_ ) return false;
  
  // Lock the src layer, so it cannot be used else where
  if ( ! ( algo->lock_for_use( algo->src_layer_ ) ) )
  {
    return false;
  }
  
  for (size_t i = 0; i < this->initalGuessSet_.size(); ++i)
  {
    LayerHandle maskLayer = LayerManager::FindLayer( this->initalGuessSet_[i], this->sandbox_ );
    if (maskLayer)
    {
      algo->initialGuessSet_.push_back(maskLayer);
    }
  }

  // Lock the mask layers, so it cannot be used else where
  for (size_t i = 0; i < algo->initialGuessSet_.size(); ++i)
  {
    if ( ! ( algo->lock_for_use( algo->initialGuessSet_[i] ) ) )
    {
      return false;
    }
  }

	algo->dst_layer_.resize( LAYER_COUNT+1 );
	
	// Create the destination layer, which will show progress
	std::vector< std::string > dst_layer_ids( LAYER_COUNT+1 );
  // hardcoded...
	for ( size_t j = 0; j < LAYER_COUNT;  j++ )
	{
		algo->create_and_lock_mask_layer_from_layer( algo->src_layer_, algo->dst_layer_[ j ] );
		dst_layer_ids.push_back( algo->dst_layer_[ j ]->get_layer_id() );
	}
	
	// Return the id of the destination layer.
	if ( algo->dst_layer_.size() > 0 )
	{
		result = Core::ActionResultHandle( new Core::ActionResult( dst_layer_ids ) );
	}
  
  // If the action is run from a script (provenance is a special case of script),
  // return a notifier that the script engine can wait on.
  if ( context->source() == Core::ActionSource::SCRIPT_E ||
       context->source() == Core::ActionSource::PROVENANCE_E )
  {
    context->report_need_resource( algo->get_notifier() );
  }

  // Build the undo-redo record
  algo->create_undo_redo_and_provenance_record( context, this->shared_from_this() );
  
  // Start the filter on a separate thread.
  Core::Runnable::Start( algo );
  
  return true;
}

void ActionReconstructionTool::Dispatch( Core::ActionContextHandle context,
                                         std::string target_layer,
                                         const std::vector< std::string >& initialGuessSet,
                                         std::string outputDir,
                                         int iterations,
                                         double xyVoxelSizeScale,
                                         double zVoxelSizeScale,
                                         callback_type callback )
{
  ActionReconstructionTool* action = new ActionReconstructionTool;
  
  action->target_layer_ = target_layer;
  action->initalGuessSet_ = initialGuessSet;
  action->outputDir_ = outputDir;
  action->iterations_ = iterations;
  action->xyVoxelSizeScale_ = xyVoxelSizeScale;
  action->zVoxelSizeScale_ = zVoxelSizeScale;
  action->callback_ = callback;

  // hack to force sandbox reset before running algorithm
  action->sandbox_ = -1;

  Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

void ActionReconstructionTool::clear_cache()
{
  std::cerr << "ActionReconstructionTool::clear_cache()" << std::endl;
  // Reset the callback so it won't be holding any handles
  this->callback_ = 0;
}

} // end namespace Seg3D
