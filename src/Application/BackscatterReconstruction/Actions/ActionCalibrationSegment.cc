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

#include <Application/BackscatterReconstruction/Actions/ActionCalibrationSegment.h>
#include <Application/BackscatterReconstruction/Algorithm/recon_api.h>
#include <Application/Filters/ITKFilter.h>

#include <Core/DataBlock/MaskDataBlockManager.h>

#include <sstream>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, CalibrationSegment )

#define LAYER_COUNT 3

namespace Seg3D
{

class CalibrationSegmentAlgo : public ITKFilter
{
  
public:
  LayerHandle src_layer_;
	std::vector<LayerHandle> dst_layer_;

public:
  // RUN:
  // Implementation of run of the Runnable base class, this function is called when the thread
  // is launched.
  
  // NOTE: The macro needs a data type to select which version to run. This needs to be
  // a member variable of the algorithm class.
  SCI_BEGIN_ITK_RUN()
  {
		// Retrieve the image as an itk image from the underlying data structure
		// NOTE: This only does wrapping and does not regenerate the data.
		Core::ITKImageDataT<float>::Handle image; 
		this->get_itk_image_from_layer<float>( this->src_layer_, image );
    
    // size set in SegmentCalibPoints
    std::vector<unsigned char> diskIds;
    
    CalibrationSegment(image->get_image(), diskIds);

    for (size_t i = 0; i < LAYER_COUNT; ++i)
    {
      MaskLayerHandle mask_layer = boost::dynamic_pointer_cast<MaskLayer>( this->dst_layer_[ i ] );
      Core::MaskDataBlockHandle maskDataBlock;
      Core::MaskDataBlockManager::Create( this->dst_layer_[ i ]->get_grid_transform(), maskDataBlock );
      Core::MaskVolumeHandle mask_volume( new Core::MaskVolume(this->dst_layer_[ i ]->get_grid_transform(), maskDataBlock ) );

      for (size_t j = 0; j < maskDataBlock->get_size(); j++ )
      {
        if ( diskIds[j] == (i+1) )
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
    return "CalibrationSegmentAlgorithm";
  }
  
  // GET_LAYER_PREFIX:
  // This function returns the name of the filter. The latter is prepended to the new layer name, 
  // when a new layer is generated. 
  virtual std::string get_layer_prefix() const
  {
    return "CalibrationSegment";
  }
};

bool ActionCalibrationSegment::validate( Core::ActionContextHandle& context )
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

bool ActionCalibrationSegment::run( Core::ActionContextHandle& context,
                                Core::ActionResultHandle& result )
{
	// Create algorithm
	boost::shared_ptr<CalibrationSegmentAlgo> algo( new CalibrationSegmentAlgo );
  
	// Copy the parameters over to the algorithm that runs the filter
	algo->set_sandbox( this->sandbox_ );

	// Find the handle to the layer
	algo->src_layer_ = LayerManager::FindLayer( this->target_layer_, this->sandbox_ );
	// Check if layer really exists
	if ( ! algo->src_layer_ ) return false;
  
	// Lock the src layer, so it cannot be used else where
	if ( ! ( algo->lock_for_use( algo->src_layer_ ) ) )
	{
		return false;
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
  
  
void ActionCalibrationSegment::Dispatch( Core::ActionContextHandle context, std::string target_layer )
{
  ActionCalibrationSegment* action = new ActionCalibrationSegment;
	
  action->target_layer_ = target_layer;

  Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}
  
} // end namespace Seg3D
