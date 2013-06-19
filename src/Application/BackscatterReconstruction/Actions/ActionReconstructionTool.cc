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
#include <Core/Utils/Log.h>

#include <Core/DataBlock/ITKDataBlock.h>
#include <Core/DataBlock/ITKImageData.h>
#include <Core/DataBlock/MaskDataBlockManager.h>

// Application Includes
#include <Application/BackscatterReconstruction/Actions/ActionReconstructionTool.h>
#include <Application/BackscatterReconstruction/Algorithm/recon_api.h>

#include <Application/ToolManager/Actions/ActionOpenTool.h>
#include <Application/ViewerManager/ViewerManager.h>
#include <Application/Layer/LayerManager.h>
#include <Application/Layer/LayerMetaData.h>

//#include <Application/Filters/ITKFilter.h>
//#include <Application/BackscatterReconstruction/Filters/ReconstructionFilter.h>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/conversion.hpp>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, ReconstructionTool )

namespace Seg3D
{

class ReconstructionToolAlgo : public ReconstructionFilter
{  

public:
  LayerHandle src_layer_;
  std::vector< LayerHandle > initialGuessSet_;

//  std::vector<LayerHandle> dst_layer_;

  std::string outputDir_;
  double xyVoxelSizeScale_;
  double zVoxelSizeScale_;
  int iterations_;

  ReconstructionToolAlgo(ReconstructionFilter::progress_callback callback, int layerCount)
    : ReconstructionFilter(callback, layerCount)
  {}
  
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
    ReconstructionFilter::FLOAT_IMAGE_TYPE::SizeType inSize = image->get_image()->GetLargestPossibleRegion().GetSize();

    DataLayerHandle data_layer = boost::dynamic_pointer_cast<DataLayer>( this->src_layer_ );
    Core::VolumeHandle srcVolumeHandle = data_layer->get_volume();
    
    ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer initialGuess;

    if (this->initialGuessSet_.size() > 0)
    {
      initialGuess = ReconstructionFilter::UCHAR_IMAGE_TYPE::New();
      ReconstructionFilter::UCHAR_IMAGE_TYPE::SizeType size;
      size.SetElement(0, inSize[0]);
      size.SetElement(1, inSize[1]);
      size.SetElement(2, inSize[2]);
      initialGuess->SetRegions(ReconstructionFilter::UCHAR_IMAGE_TYPE::RegionType(size));
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
        }
      }
    }

    float voxelSizeCM[3] = { this->xyVoxelSizeScale_, this->xyVoxelSizeScale_, this->zVoxelSizeScale_ };

    Layer::filter_key_type key = this->get_key();
    SandboxID sandbox = this->get_sandbox();
    this->start_progress();

    int test_iteration_count = 0;
    ReconstructionStart(image->get_image(),
                        initialGuess,
                        algorithm_config_file.string().c_str(),
                        voxelSizeCM,
                        /*iterations_*/ test_iteration_count,
                        get_recon_volume());

    this->stop_progress();
    
//    ReconstructionFilter::UCHAR_IMAGE_TYPE::SizeType outSize = finalReconVolume->GetLargestPossibleRegion().GetSize();
//std::cerr << "Output ITK image size: " << outSize[0] << ", " << outSize[1] << ", " << outSize[2] << std::endl;
//    const size_t SIZE = outSize[0] * outSize[1] * outSize[2];
//    if (SIZE == 0)
//    {
//      CORE_LOG_WARNING("Reconstruction mask size is 0.");
//      return;
//    }
//    Core::DataBlockHandle finalReconDataBlock = Core::ITKDataBlock::New(finalReconVolume.GetPointer());
//    Core::ITKUCharImageDataHandle finalReconImage = Core::ITKUCharImageDataHandle(
//      new Core::ITKUCharImageData(finalReconDataBlock) );
//
//    MaskLayerHandle maskLayer = boost::dynamic_pointer_cast<MaskLayer>( this->dst_layer_[ 0 ] );
//    Core::DataBlockHandle dataBlock = Core::ITKDataBlock::New(finalReconVolume.GetPointer());
//    Core::ITKUCharImageDataHandle outImage = Core::ITKUCharImageDataHandle(
//      new Core::ITKUCharImageData(dataBlock) );
//    
//    this->dst_layer_[0]->set_grid_transform(outImage->get_grid_transform(), false);
//    
//    Core::MaskDataBlockHandle maskDataBlock;
//    if (!( Core::MaskDataBlockManager::Convert( dataBlock, outImage->get_grid_transform(), maskDataBlock ) ) )
//    {
//      CORE_LOG_WARNING("Could not allocate enough memory for temporary mask.");
//    }
//    Core::MaskVolumeHandle maskVolume( new Core::MaskVolume(outImage->get_grid_transform(), maskDataBlock ) );
//
//    for (size_t j = 0; j < maskDataBlock->get_size(); j++ )
//    {
//      // test - insert all
//      if ( finalReconVolume->GetBufferPointer()[j] > 0 )
//      {
//        maskDataBlock->set_mask_at(j);
//      }
//    }
//    
//    if (! this->dispatch_insert_mask_volume_into_layer( this->dst_layer_[ 0 ], maskVolume ) )
//    {
//      std::ostringstream oss;
//      oss << "Could not insert mask volume " << this->dst_layer_[ 0 ]->get_layer_id() << " into layer";
//      this->report_error(oss.str());
//    }
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
  // TODO: more systematic way to know number of masks in output?
  const int LAYER_COUNT = 2;
  // Create algorithm
  boost::shared_ptr<ReconstructionToolAlgo> algo( new ReconstructionToolAlgo(this->callback_, LAYER_COUNT) );

  // Copy the parameters over to the algorithm that runs the filter
  algo->set_sandbox( this->sandbox_ );

  // Set up parameters
  algo->outputDir_ = this->outputDir_;
  algo->xyVoxelSizeScale_ = this->xyVoxelSizeScale_;
  algo->zVoxelSizeScale_ = this->zVoxelSizeScale_;
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
                                         ReconstructionFilter::progress_callback callback )
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

void ActionReconstructionTool::Abort()
{
  ReconstructionAbort();
}

//void ActionReconstructionTool::Test(MaskLayerHandle layer, Core::MaskVolumeHandle mask, ProvenanceID prov_id, Layer::filter_key_type key, SandboxID sandbox)
//{
//  LayerManager::DispatchInsertMaskVolumeIntoLayer( layer, mask, prov_id, key, sandbox );
//}


void ActionReconstructionTool::clear_cache()
{
  std::cerr << "ActionReconstructionTool::clear_cache()" << std::endl;
  // Reset the callback so it won't be holding any handles
  this->callback_ = 0;
}

} // end namespace Seg3D
