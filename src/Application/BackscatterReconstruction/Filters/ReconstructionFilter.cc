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

class ReconstructionFilterProgress : public Core::Runnable
{
public:
  LayerHandle tmpLayer_;
 
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
    tmpLayer_.reset();
  }
  
  void stop() { stop_ = true; }
  
  void setTemporaryLayer(const LayerHandle& layerHandle) { this->tmpLayer_ = layerHandle; }
  LayerHandle getTemporaryLayer() { return this->tmpLayer_; }
  
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
std::cerr << "Temp ITK image size: " << sizes[0] << ", " << sizes[1] << ", " << sizes[2] << std::endl;
      if (SIZE > 0)
      {
        Core::DataBlockHandle tmp_dataBlock = Core::ITKDataBlock::New(reconVolume.GetPointer());
        Core::ITKUCharImageDataHandle tmp_image = Core::ITKUCharImageDataHandle(
          new Core::ITKUCharImageData(tmp_dataBlock) );
        
        if (! tmpLayer_ && ! layerCreateFlag_)
        {
          LayerMetaData metadata;
          Core::Application::PostEvent( boost::bind(
            &LayerManager::CreateAndLockMaskLayer, tmp_image->get_grid_transform(),
            ReconstructionFilter::TMP_LAYER_NAME, tmpLayer_,
            metadata, this->filter_->get_key(), this->filter_->get_sandbox() ) );
          // TODO: check for completions????
          layerCreateFlag_ = true;
          
          continue;
        }
        
        if (tmpLayer_)
        {
          std::cerr << "Temporary Reconstruction layer created" << std::endl;
          Core::MaskDataBlockHandle maskDataBlock;
          if (!( Core::MaskDataBlockManager::Convert( tmp_dataBlock, tmp_image->get_grid_transform(), maskDataBlock ) ) )
          {
            CORE_LOG_WARNING("Could not allocate enough memory for temporary mask.");
          }
          else
          {
std::cerr << "Attempt to insert data into layer..." << std::endl;
            // move to filter, add filter support...
            Core::MaskVolumeHandle maskVolume( new Core::MaskVolume(tmp_image->get_grid_transform(), maskDataBlock ) );
            MaskLayerHandle tmpMaskLayer = boost::dynamic_pointer_cast<MaskLayer>( tmpLayer_ );
            // not bothering with provenance here...
            Core::Application::PostEvent( boost::bind(
              &LayerManager::DispatchInsertMaskVolumeIntoLayer, tmpMaskLayer, maskVolume, -1,
                this->filter_->get_key(), this->filter_->get_sandbox() ) );
          }
        }

//        if ( tmpLayer_ && tmpLayer_->locked_state_->get() )
//        {
//std::cerr << "Attempt unlock layer" << std::endl;
//          Core::Application::PostEvent( boost::bind(
//            &LayerManager::DispatchUnlockLayer, tmpLayer_, key_, sandbox_ ) );
//        }
//        else
//        {
//          std::cerr << "Layer unavailable..." << std::endl;
//          if (! tmpLayer_ )
//          {
//            std::cerr << "Layer is NULL" << std::endl;
//          }
//          else if ( tmpLayer_->locked_state_->get() )
//          {
//            std::cerr << "Layer is locked" << std::endl;
//          }
//        }
//        
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
      layerCount_(layerCount)
  {
    dstLayers_.resize(this->layerCount_+1);
  }

  ~ReconstructionFilterPrivate()
  {
  }
  
	void handle_layer_group_insert( LayerHandle layerHandle, bool newGroup );
  void finalize();

  LayerHandle tmpLayer_;
  std::vector<LayerHandle> dstLayers_;

  ReconstructionFilterProgressHandle progress_;
  ReconstructionFilterHandle filter_;
  ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer finalReconVolume_;
  int layerCount_;
};

void ReconstructionFilterPrivate::handle_layer_group_insert( LayerHandle layerHandle, bool newGroup )
{
  lock_type lock( this->get_mutex() );

  std::cerr << "Detected " << layerHandle->get_layer_id() << ", " << layerHandle->get_layer_name() << std::endl;
  if (layerHandle->get_layer_name() == ReconstructionFilter::TMP_LAYER_NAME)
  {
    this->progress_->tmpLayer_ = layerHandle;
  }

  if (newGroup)
  {
    for (int i = 0; i < this->layerCount_; ++i)
    {
      LayerMetaData metadata;
      Core::Application::PostEvent( boost::bind(
        &LayerManager::CreateAndLockMaskLayer, layerHandle->get_grid_transform(), "test", 
        this->dstLayers_[i], metadata, this->filter_->get_key(), this->filter_->get_sandbox() ) );
    }
  }
}

void ReconstructionFilterPrivate::finalize()
{
std::cerr << "ReconstructionFilterPrivate::finalize() begin" << std::endl;
  
  ReconstructionFilter::UCHAR_IMAGE_TYPE::SizeType outSize = finalReconVolume_->GetLargestPossibleRegion().GetSize();
  std::cerr << "Output ITK image size: " << outSize[0] << ", " << outSize[1] << ", " << outSize[2] << std::endl;
  const size_t SIZE = outSize[0] * outSize[1] * outSize[2];

  for (size_t i = 0; i < SIZE; ++i)
  {
    if ( finalReconVolume_->GetBufferPointer()[i] > 0 )
    {
      fprintf(stderr, "recon volume mask at %lu = %hhu\n", i, finalReconVolume_->GetBufferPointer()[i]);
    }
  }

  if ( tmpLayer_ )
  {
    Core::Application::PostEvent( boost::bind(
      &LayerManager::DispatchDeleteLayer, tmpLayer_, this->filter_->get_key(), this->filter_->get_sandbox() ) );
    tmpLayer_.reset();
  }    
  this->disconnect_all();
  std::cerr << "ReconstructionFilterPrivate::finalize() end" << std::endl;
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
if ( !( Core::Application::IsApplicationThread() ) )
{
std::cerr << "NOT on Application thread" << std::endl;
}

std::cerr << "ReconstructionFilter::~ReconstructionFilter() begin" << std::endl;
  ReconstructionFilterPrivate::lock_type lock( this->private_->get_mutex() );
  this->private_->finalize();
std::cerr << "ReconstructionFilter::~ReconstructionFilter() end" << std::endl;
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
  this->private_->progress_->filter_.reset(this);
  Core::Runnable::Start( this->private_->progress_ );
}

void ReconstructionFilter::stop_progress()
{
  this->private_->progress_->stop();
}

ReconstructionFilter::UCHAR_IMAGE_TYPE::Pointer ReconstructionFilter::get_recon_volume()
{
  return this->private_->finalReconVolume_;
}

//void ReconstructionFilter::limit_number_of_itk_threads_internal( itk::ProcessObject::Pointer filter )
//{
//  // Assume we will have a minimum of 2 threads. As we subtract one this will ensure that
//  // there is at least one thread doing the computation.
//  unsigned int max_threads = boost::thread::hardware_concurrenReconstructionFiltercy();
//  if ( max_threads < 2 ) max_threads = 2;
//  
//  filter->GetMultiThreader()->SetNumberOfThreads( max_threads - 1 );
//}

} // end namespace Core
