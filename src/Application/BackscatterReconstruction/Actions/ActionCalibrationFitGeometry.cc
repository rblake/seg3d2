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

#include <Application/BackscatterReconstruction/Actions/ActionCalibrationFitGeometry.h>
#include <Application/BackscatterReconstruction/Algorithm/recon_api.h>
#include <Application/Filters/ITKFilter.h>

#include <sstream>

#define CALIBRATION_SET_SIZE 3

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, CalibrationFitGeometry )


namespace Seg3D
{

class CalibrationFitGeometryAlgo : public ITKFilter
{
  
public:
  LayerHandle src_layer_;
	//LayerHandle dst_layer_;

  std::vector< std::string > calibrationSet_;

public:
  // RUN:
  // Implementation of run of the Runnable base class, this function is called when the thread
  // is launched.
  
  // NOTE: The macro needs a data type to select which version to run. This needs to be
  // a member variable of the algorithm class.
  SCI_BEGIN_ITK_RUN()
  {
    if (this->calibrationSet_.size() < CALIBRATION_SET_SIZE)
    {
      this->report_error("Calibration ids not set");
      return;
    }
    boost::filesystem::path algorithm_work_dir;
    boost::filesystem::path algorithm_config_file;    
    boost::filesystem::path algorithm_geometry_file;
    if(! Core::Application::Instance()->get_algorithm_config(algorithm_work_dir,
           algorithm_config_file, algorithm_geometry_file) )
	{
      this->report_error("Failed to get algorithm config information.");
      return;
	}
    
    // Retrieve the image as an itk image from the underlying data structure
    // NOTE: This only does wrapping and does not regenerate the data.
    Core::ITKImageDataT<float>::Handle image; 
    this->get_itk_image_from_layer<float>( this->src_layer_, image );
    
    ITKFilter::FLOAT_IMAGE_TYPE::SizeType inSize = image->get_image()->GetLargestPossibleRegion().GetSize();
    
    DataLayerHandle data_layer = boost::dynamic_pointer_cast<DataLayer>( this->src_layer_ );
    Core::VolumeHandle srcVolumeHandle = data_layer->get_volume();
    if ( inSize[0] * inSize[1] * inSize[2] != srcVolumeHandle->get_size())
    {
      this->report_error("ITK region size does not match data volume size");
      return;
    }

    std::vector<unsigned char> maskv(inSize[0] * inSize[1] * inSize[2]);
    
    for (size_t i = 0; i < this->calibrationSet_.size(); ++i)
    {
      LayerHandle layer = LayerManager::FindLayer( this->calibrationSet_[i], this->get_sandbox() );
      if (! layer)
      {
        std::ostringstream oss;
        oss << "Could not find layer " << this->calibrationSet_[i];
        this->report_error(oss.str());
        return;
      }

      // TODO: lock layer?

      Core::MaskLayerHandle mask_layer = boost::dynamic_pointer_cast<MaskLayer>( layer );
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
          // set disk IDs from 1 to 3
          maskv[j] = (i+1);
        }
        else
        {
          maskv[j] = 0;
        }
      }
    }

    CalibrationFitGeometry(image->get_image(),
                           maskv,
                           algorithm_work_dir.string().c_str(),
                           algorithm_config_file.string().c_str(),
                           algorithm_geometry_file.string().c_str());

    // TODO: output besides new geometry file?

  }
  SCI_END_ITK_RUN()
  
  // GET_FITLER_NAME:
  // The name of the filter, this information is used for generating new layer labels.
  virtual std::string get_filter_name() const
  {
    return "CalibrationFitGeometryAlgorithm";
  }
  
  // GET_LAYER_PREFIX:
  // This function returns the name of the filter. The latter is prepended to the new layer name, 
  // when a new layer is generated. 
  virtual std::string get_layer_prefix() const
  {
    return "CalibrationFitGeometry";
  }
};

bool ActionCalibrationFitGeometry::validate( Core::ActionContextHandle& context )
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

bool ActionCalibrationFitGeometry::run( Core::ActionContextHandle& context,
                                   Core::ActionResultHandle& result )
{  
  // Create algorithm
  boost::shared_ptr<CalibrationFitGeometryAlgo> algo( new CalibrationFitGeometryAlgo );
  
  // Copy the parameters over to the algorithm that runs the filter
  //algo->set_sandbox( this->sandbox_ );
  algo->calibrationSet_ = this->calibrationSet_;

  // Find the handle to the layer
  algo->src_layer_ = LayerManager::FindLayer( this->target_layer_, this->sandbox_ );
  // Check if layer really exists
  if ( ! algo->src_layer_ ) return false;
  
  // Lock the src layer, so it cannot be used else where
  if ( ! ( algo->lock_for_use( algo->src_layer_ ) ) )
  {
    return false;
  }

//	// Create the destination layer, which will show progress
//	if ( !( algo->create_and_lock_mask_layer_from_layer( algo->src_layer_, algo->dst_layer_ ) ) )
//	{
//		return false;
//	}
//  
//	// Return the id of the destination layer.
//	result = Core::ActionResultHandle( new Core::ActionResult( algo->dst_layer_->get_layer_id() ) );  
  
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


void ActionCalibrationFitGeometry::Dispatch( Core::ActionContextHandle context,
                                             const std::string& target_layer,
                                             const std::vector< std::string >& calibrationSet )
{
  ActionCalibrationFitGeometry* action = new ActionCalibrationFitGeometry;
  action->target_layer_ = target_layer;
  action->calibrationSet_ = calibrationSet;
  // hack to force sandbox reset before running algorithm
  action->sandbox_ = -1;
  
  Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}
  
} // end namespace Seg3D
