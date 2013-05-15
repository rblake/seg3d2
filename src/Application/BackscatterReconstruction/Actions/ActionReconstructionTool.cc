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

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, ReconstructionTool )

namespace Seg3D
{

class ReconstructionToolAlgo : public ITKFilter
{  

public:
  LayerHandle src_layer_;
  std::vector<LayerHandle> dst_layer_;

  std::vector< std::string > initialGuessSet_;
  std::string outputDir_;
  double measurementScale_;
  int iterations_;
  
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
    
		Core::ITKImageDataT<float>::Handle image; 
		this->get_itk_image_from_layer<float>( this->src_layer_, image );
    
    // null for now...
    ITKFilter::UCHAR_IMAGE_TYPE::Pointer initialGuess = ITKFilter::UCHAR_IMAGE_TYPE::New();  

    // TODO: not implemented, so returns immediately
    ReconstructionStart(image->get_image(),
                        initialGuess,
                        algorithm_config_file.string().c_str(),
                        iterations_);
    // TODO: will have to bind tool signal to API function
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
  this->add_parameter( this->iterations_ );
  this->add_parameter( this->outputDir_ );
  this->add_parameter( this->measurementScale_ );
  this->add_parameter( this->initalGuessSet_ );
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
  algo->initialGuessSet_ = this->initalGuessSet_;
  algo->measurementScale_ = this->measurementScale_;
  algo->iterations_ = this->iterations_;
  
	// Find the handle to the layer
	algo->src_layer_ = LayerManager::FindLayer( this->target_layer_, this->sandbox_ );
	// Check if layer really exists
	if ( ! algo->src_layer_ ) return false;
  
	// Lock the src layer, so it cannot be used else where
	if ( ! ( algo->lock_for_use( algo->src_layer_ ) ) )
	{
		return false;
	}
  
  // TODO: set up destination layers once output is known
  
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
                                         std::string outputDir,
                                         const std::vector< std::string >& initialGuessSet,
                                         int iterations,
                                         double measurementScale )
{
  ActionReconstructionTool* action = new ActionReconstructionTool;
  
  action->target_layer_ = target_layer;
  action->outputDir_ = outputDir;
  action->initalGuessSet_ = initialGuessSet;
  action->iterations_ = iterations;
  action->measurementScale_ = measurementScale;
  // hack to force sandbox reset before running algorithm
  action->sandbox_ = -1;

  Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}
  
} // end namespace Seg3D
