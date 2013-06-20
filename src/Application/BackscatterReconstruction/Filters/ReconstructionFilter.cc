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

// ITK includes 
#include <itkCommand.h>

// Core includes
#include <Core/Utils/Exception.h>
#include <Core/DataBlock/StdDataBlock.h>
#include <Core/DataBlock/ITKDataBlock.h>
#include <Core/DataBlock/ITKImageData.h>
#include <Core/DataBlock/MaskDataBlock.h>
#include <Core/DataBlock/MaskDataBlockManager.h>

// Application includes
#include <Application/Layer/Layer.h>
#include <Application/Filters/LayerFilterLock.h>
#include <Application/BackscatterReconstruction/Filters/ReconstructionFilter.h>
#include <Application/BackscatterReconstruction/Algorithm/recon_api.h>

#include <Core/Utils/Log.h>


namespace Seg3D
{
  
const std::string ReconstructionFilter::TMP_LAYER_NAME("TemporaryReconstruction");
const std::string ReconstructionFilter::TMP_LAYER_META_INFO("TemporaryReconstructionMaskLayer");
const std::string ReconstructionFilter::DEST_LAYER_META_INFO("ReconstructionOutputMaskLayer");

class ReconstructionFilterProgress : public Core::Runnable
{
public:
//  LayerHandle tmpLayer_;
 
  long msTimeInterval_;
  bool stop_;
  bool layerCreateFlag_;

  ReconstructionFilterHandle filter_;
  ReconstructionFilter::progress_callback callback_;
  
  ReconstructionFilterProgress(ReconstructionFilter::progress_callback callback,
                               long msTimeInterval=500)
  : callback_(callback),
    msTimeInterval_(msTimeInterval),
    stop_(false),
    layerCreateFlag_(false) {}
  
  ~ReconstructionFilterProgress()
  {
  }
  
  void stop() { stop_ = true; }
  
  virtual void run()
  {
    // TODO: temporary? See TODO below...
    double scale = 100;
    while (! stop_ )
    {
      ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume = ReconstructionFilter::UCHAR_IMAGE_TYPE::New();
      double progress = ReconstructionGetMaterialVolume(reconVolume);
      if ( this->callback_ != 0 )
      {
        // TODO: using scale temporarily because progress values are < 0.01
        // before reconstruction either exits or crashes
        Core::Application::PostEvent( boost::bind( this->callback_, progress*scale, 0, 1.0 ) );
      }
      // test
      ReconstructionFilter::UCHAR_IMAGE_TYPE::SizeType sizes = reconVolume->GetLargestPossibleRegion().GetSize();
      const size_t SIZE = sizes[0] * sizes[1] * sizes[2];
      if (SIZE > 0)
      {
        if (! layerCreateFlag_)
        {
          this->filter_->create_and_lock_tmp_mask_layer(reconVolume);
          layerCreateFlag_ = true;
          
          continue;
        }
        else
        {
          this->filter_->update_tmp_mask_layer(reconVolume);
        }
      }

      boost::this_thread::sleep(boost::posix_time::milliseconds(msTimeInterval_));
    }
    // make sure progress bar completes
    Core::Application::PostEvent( boost::bind( this->callback_, 100.0, 0, 1.0 ) );
  }
};


class ReconstructionFilterPrivate : public Core::Lockable, public Core::ConnectionHandler
{
public:
  ReconstructionFilterPrivate(ReconstructionFilter::progress_callback callback, int layerCount)
    : progress_( new ReconstructionFilterProgress(callback) ),
      finalReconVolume_( ReconstructionFilter::UCHAR_IMAGE_TYPE::New() ),
      layerCount_(layerCount),
      layerName_("reconstruction")
  {
    dstLayers_.resize(this->layerCount_+1);
  }

  ~ReconstructionFilterPrivate()
  {
    tmpLayer_.reset();
  }
  
	void handle_layer_group_insert( LayerHandle layerHandle, bool newGroup );
  void create_tmp_layer(ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume);
  void update_tmp_layer(ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume);

  void finalize();

  LayerHandle tmpLayer_;
  std::vector<LayerHandle> dstLayers_;

  ReconstructionFilterProgressHandle progress_;
  ReconstructionFilterHandle filter_;
  ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer finalReconVolume_;
  int layerCount_;
  std::string layerName_;
};

void ReconstructionFilterPrivate::update_tmp_layer(ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume)
{
  if (tmpLayer_)
  {
    Core::DataBlockHandle tmp_dataBlock = Core::ITKDataBlock::New(reconVolume.GetPointer());
    Core::ITKUCharImageDataHandle tmp_image = Core::ITKUCharImageDataHandle(
      new Core::ITKUCharImageData(tmp_dataBlock) );

    Core::MaskDataBlockHandle maskDataBlock;
    if (!( Core::MaskDataBlockManager::Convert( tmp_dataBlock, tmp_image->get_grid_transform(), maskDataBlock ) ) )
    {
      CORE_LOG_WARNING("Could not allocate enough memory for temporary mask.");
    }
    else
    {
      std::ostringstream oss;
      oss << "Attempt to insert progress data into temporary layer...";
      CORE_LOG_MESSAGE(oss.str());

      // move to filter, add filter support...
      Core::MaskVolumeHandle maskVolume( new Core::MaskVolume(tmp_image->get_grid_transform(), maskDataBlock ) );
      MaskLayerHandle tmpMaskLayer = boost::dynamic_pointer_cast<MaskLayer>( tmpLayer_ );
      // not bothering with provenance here...
      Core::Application::PostEvent( boost::bind(
        &LayerManager::DispatchInsertMaskVolumeIntoLayer, tmpMaskLayer, maskVolume, -1,
        this->filter_->get_key(), this->filter_->get_sandbox() ) );
    }
  }
}

void ReconstructionFilterPrivate::create_tmp_layer(ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume)
{
  if (! tmpLayer_ )
  {
    Core::DataBlockHandle tmp_dataBlock = Core::ITKDataBlock::New(reconVolume.GetPointer());
    Core::ITKUCharImageDataHandle tmp_image = Core::ITKUCharImageDataHandle(
      new Core::ITKUCharImageData(tmp_dataBlock) );

    LayerMetaData metadata;
    metadata.meta_data_info_ = ReconstructionFilter::TMP_LAYER_META_INFO;
    metadata.meta_data_ = ReconstructionFilter::TMP_LAYER_NAME;

    Core::Application::PostEvent( boost::bind(
      &LayerManager::CreateAndLockMaskLayer, tmp_image->get_grid_transform(),
      ReconstructionFilter::TMP_LAYER_NAME, tmpLayer_,
      metadata, this->filter_->get_key(), this->filter_->get_sandbox() ) );
  }
}

void ReconstructionFilterPrivate::handle_layer_group_insert( LayerHandle layerHandle, bool newGroup )
{
  lock_type lock( this->get_mutex() );

  std::ostringstream oss1;
  oss1 << "Detected " << layerHandle->get_layer_id() << ", " << layerHandle->get_layer_name() << " in handler";
  CORE_LOG_MESSAGE(oss1.str());

  if (layerHandle->get_layer_name() == ReconstructionFilter::TMP_LAYER_NAME &&
      layerHandle->get_meta_data().meta_data_info_ == ReconstructionFilter::TMP_LAYER_META_INFO)
  {
    this->tmpLayer_ = layerHandle;

    if (newGroup)
    {
      for (int i = 0; i < this->layerCount_; ++i)
      {
        LayerMetaData metadata;
        metadata.meta_data_info_ = ReconstructionFilter::DEST_LAYER_META_INFO;
        std::ostringstream oss;
        oss << i;
        metadata.meta_data_ = oss.str();

        Core::Application::PostEvent( boost::bind(
          &LayerManager::CreateAndLockMaskLayer, layerHandle->get_grid_transform(), this->layerName_, 
          this->dstLayers_[i], metadata, this->filter_->get_key(), this->filter_->get_sandbox() ) );
      }
    }
  }
  else if (layerHandle->get_meta_data().meta_data_info_ == ReconstructionFilter::DEST_LAYER_META_INFO)
  {
    std::istringstream iss(layerHandle->get_meta_data().meta_data_);
    int index;
    iss >> index;
    std::ostringstream oss2;
    oss2 << "Processing destination layer " << index << " from handler";
    CORE_LOG_MESSAGE(oss2.str());

    if (index >= 0 && index < this->layerCount_)
    {
      this->dstLayers_[index] = layerHandle;
    }
  }
}

void ReconstructionFilterPrivate::finalize()
{
  if (this->tmpLayer_)
  {
//    Core::Application::PostEvent( boost::bind(
//      &LayerManager::DispatchUnlockOrDeleteLayer, tmpLayer_, this->filter_->get_key(), this->filter_->get_sandbox() ) );
    Core::Application::PostEvent( boost::bind(
      &LayerManager::DispatchDeleteLayer, tmpLayer_, this->filter_->get_key(), this->filter_->get_sandbox() ) );
  }
  
  ReconstructionFilter::UCHAR_IMAGE_TYPE::SizeType outSize = finalReconVolume_->GetLargestPossibleRegion().GetSize();

  const size_t SIZE = outSize[0] * outSize[1] * outSize[2];
  if (SIZE > 0)
  {
    Core::DataBlockHandle finalReconDataBlock = Core::ITKDataBlock::New(this->finalReconVolume_.GetPointer());
    Core::ITKUCharImageDataHandle finalReconImage = Core::ITKUCharImageDataHandle(
      new Core::ITKUCharImageData(finalReconDataBlock) );

    for (size_t i = 0; i < this->layerCount_; ++i)
    {
      MaskLayerHandle maskLayer = boost::dynamic_pointer_cast<MaskLayer>( this->dstLayers_[ i ] );

      Core::DataBlockHandle dataBlock = Core::ITKDataBlock::New(this->finalReconVolume_.GetPointer());
      Core::ITKUCharImageDataHandle outImage = Core::ITKUCharImageDataHandle(
        new Core::ITKUCharImageData(dataBlock) );
      
      Core::MaskDataBlockHandle maskDataBlock;
      if (! ( Core::MaskDataBlockManager::Convert( dataBlock, outImage->get_grid_transform(), maskDataBlock ) ) )
      {
        CORE_LOG_WARNING("Could not allocate enough memory for temporary mask.");
      }

      Core::MaskVolumeHandle maskVolume( new Core::MaskVolume(outImage->get_grid_transform(), maskDataBlock ) );
      for (size_t j = 0; j < maskDataBlock->get_size(); ++j )
      {
        if ( this->finalReconVolume_->GetBufferPointer()[j] == (i+1) )
        {
          maskDataBlock->set_mask_at(j);
        }
      }

      // not bothering with provenance here...
      Core::Application::PostEvent( boost::bind(
        &LayerManager::DispatchInsertMaskVolumeIntoLayer, maskLayer, maskVolume, -1,
        this->filter_->get_key(), this->filter_->get_sandbox() ) );

      Core::Application::PostEvent( boost::bind(
        &LayerManager::DispatchUnlockLayer, maskLayer,
        this->filter_->get_key(), this->filter_->get_sandbox() ) );
    }
  }
  else
  {
    CORE_LOG_WARNING("Reconstruction mask size is 0.");
  }

  this->disconnect_all();
}

ReconstructionFilter::ReconstructionFilter(progress_callback callback, int layerCount) :
  LayerFilter(),
  private_( new ReconstructionFilterPrivate(callback, layerCount) )
{
  this->private_->filter_.reset(this);
  this->private_->add_connection( LayerManager::Instance()->layer_inserted_signal_.connect( 
    boost::bind( &ReconstructionFilterPrivate::handle_layer_group_insert, this->private_, _1, _2 ) ) );
}

ReconstructionFilter::~ReconstructionFilter()
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->finalize();
}

//
//void ReconstructionFilter::forward_abort_to_filter_internal( itk::ProcessObject::Pointer filter, 
//                                                 LayerHandle layer )
//{		
//  // TODO: Most of this can be handled by the layer above --JGS
//  
//  // Setup forwarding of the abort signal to the itk filter
//  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
//  this->private_->disconnect_all();
//  
//  // NOTE: The following logic is already done by LayerFilter.
//  //this->private_->add_connection( layer->abort_signal_.connect( boost::bind(
//  //	&ReconstructionFilter::raise_abort, this ) ) );
//  //this->private_->add_connection( layer->stop_signal_.connect( boost::bind(
//  //	&ReconstructionFilter::raise_stop, this ) ) );
//}

void ReconstructionFilter::handle_abort()
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
//  if ( this->private_->filter_.GetPointer() )
//  {
//    this->private_->filter_->SetAbortGenerateData( true );
//  }
}

void ReconstructionFilter::handle_stop()
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
//  if ( this->private_->filter_.GetPointer() )
//  {
//    this->private_->filter_->SetAbortGenerateData( true );
//  }
}

void ReconstructionFilter::start_progress()
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->progress_->filter_.reset(this);
  Core::Runnable::Start( this->private_->progress_ );
}

void ReconstructionFilter::stop_progress()
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->progress_->stop();
}

ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer ReconstructionFilter::get_recon_volume()
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  return this->private_->finalReconVolume_;
}

void ReconstructionFilter::set_layer_name(const std::string& name)
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->layerName_ = name;
}

void ReconstructionFilter::create_and_lock_tmp_mask_layer(ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume)
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->create_tmp_layer(reconVolume);
}

void ReconstructionFilter::update_tmp_mask_layer(ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer reconVolume)
{
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->update_tmp_layer(reconVolume);
}

} // end namespace Core
