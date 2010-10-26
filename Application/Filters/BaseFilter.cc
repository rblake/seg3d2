/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2009 Scientific Computing and Imaging Institute,
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
 
// STL includes
#include <vector> 
 
// Application includes
#include <Application/LayerManager/LayerManager.h>
#include <Application/Filters/BaseFilter.h>
#include <Application/Filters/BaseFilterLock.h>
#include <Application/StatusBar/StatusBar.h>
 
namespace Seg3D
{

class BaseFilterPrivate : public Core::ConnectionHandler
{
public:
	// Keep track of errors
	std::string error_;
	
	// Keep track of which layers were locked.
	std::vector<LayerHandle> locked_layers_;
	
	// Keep track of which layers were created.
	std::vector<LayerHandle> created_layers_;
	
	// Keep track of abort status
	bool abort_;
	
	// Mutex protecting abort status
	boost::mutex abort_mutex_;

	// Key used for this filter
	Layer::data_processing_key_type key_;

	// Function for handling abort
	void handle_abort();
};

void BaseFilterPrivate::handle_abort()
{
	boost::mutex::scoped_lock lock( abort_mutex_ );
	abort_ = true;
}

BaseFilter::BaseFilter() :
	private_( new BaseFilterPrivate )
{
	this->private_->abort_ = false;
	this->private_->key_ = Layer::GenerateDataProcessingKey();
}

BaseFilter::~BaseFilter()
{
	this->private_->disconnect_all();

	for ( size_t j = 0; j < this->private_->locked_layers_.size(); j++ )
	{
		LayerManager::DispatchUnlockLayer( this->private_->locked_layers_[ j ],
			this->private_->key_ );
	}

	boost::mutex::scoped_lock lock( this->private_->abort_mutex_ );
	if ( this->private_->abort_ )
	{
		for ( size_t j = 0; j < this->private_->created_layers_.size(); j++ )
		{
			LayerManager::DispatchDeleteLayer( this->private_->created_layers_[ j ],
				this->private_->key_ );
		}
	}
	else
	{
		for ( size_t j = 0; j < this->private_->created_layers_.size(); j++ )
		{
			LayerManager::DispatchUnlockOrDeleteLayer( this->private_->created_layers_[ j ],
				this->private_->key_ );
		}
	}
	
	
	if ( this->private_->error_.size() )
	{
		StatusBar::SetMessage( Core::LogMessageType::ERROR_E, this->private_->error_ );	
	}
}

void BaseFilter::raise_abort()
{
	boost::mutex::scoped_lock lock( this->private_->abort_mutex_ );
	this->private_->abort_ = true;
	this->report_error( "Processing was aborted." );
}

bool BaseFilter::check_abort()
{
	boost::mutex::scoped_lock lock( this->private_->abort_mutex_ );
	return this->private_->abort_;
}

void BaseFilter::connect_abort( const  LayerHandle& layer )
{
	boost::mutex::scoped_lock lock( this->private_->abort_mutex_ );
	this->private_->add_connection( layer->abort_signal_.connect( boost::bind(
		&BaseFilterPrivate::handle_abort, this->private_ ) ) );
}

void BaseFilter::report_error( const std::string& error )
{
	this->private_->error_ = this->get_filter_name() +": " + error;
}

bool BaseFilter::find_layer( const std::string& layer_id, LayerHandle& layer )
{
	layer = LayerManager::Instance()->get_layer_by_id( layer_id );
	if ( layer )
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool BaseFilter::lock_for_use( LayerHandle layer )
{
	if ( !( LayerManager::LockForUse( layer, this->private_->key_ ) ) ) 
	{
		this->report_error( "Could not lock '" + layer->get_layer_name() + "'." );
		return false;
	}
	this->private_->locked_layers_.push_back( layer );
	return true;
}

bool BaseFilter::lock_for_processing( LayerHandle layer )
{
	if ( !( LayerManager::LockForProcessing( layer, this->private_->key_ ) ) )
	{
		this->report_error( "Could not lock '" + layer->get_layer_name() + "'." );
		return false;
	}
	this->private_->locked_layers_.push_back( layer );
	return true;
}

bool BaseFilter::create_and_lock_data_layer_from_layer( LayerHandle src_layer, 
	LayerHandle& dst_layer )
{
	// Generate a new name for the filter
	std::string name = this->get_layer_prefix() + "_" + src_layer->get_layer_name();

	// Create the layer in creating mode
	if ( !( LayerManager::CreateAndLockDataLayer( src_layer->get_grid_transform(),
		name, dst_layer, this->private_->key_ ) ) )
	{
		dst_layer.reset();
		this->report_error( "Could not allocate enough memory." );
		return false;
	}
	
	// Record that the layer is locked
	this->private_->created_layers_.push_back( dst_layer );

	// Success
	return true;
}

bool BaseFilter::create_and_lock_data_layer( const Core::GridTransform& grid_trans, 
											LayerHandle src_layer, LayerHandle& dst_layer )
{
	// Generate a new name for the filter
	std::string name = this->get_layer_prefix() + "_" + src_layer->get_layer_name();

	// Create the layer in creating mode
	if ( !( LayerManager::CreateAndLockDataLayer( grid_trans, name, dst_layer,
		this->private_->key_ ) ) )
	{
		dst_layer.reset();
		this->report_error( "Could not allocate enough memory." ); 
		return false;
	}

	// Record that the layer is locked
	this->private_->created_layers_.push_back( dst_layer );

	// Success
	return true;
}

bool BaseFilter::create_and_lock_mask_layer_from_layer( LayerHandle src_layer, LayerHandle& dst_layer )
{
	// Generate a new name for the filter
	std::string name = this->get_layer_prefix() + "_" + src_layer->get_layer_name();

	// Create the layer in creating mode
	if ( !( LayerManager::CreateAndLockMaskLayer( src_layer->get_grid_transform(),
		name, dst_layer, this->private_->key_ ) ) )
	{
		dst_layer.reset();
		this->report_error( "Could not allocate enough memory." );
		return false;
	}
	
	// Record that the layer is locked
	this->private_->created_layers_.push_back( dst_layer );

	// Success
	return true;
}

bool BaseFilter::create_and_lock_mask_layer( const Core::GridTransform& grid_trans, 
											LayerHandle src_layer, LayerHandle& dst_layer )
{
	// Generate a new name for the filter
	std::string name = this->get_layer_prefix() + "_" + src_layer->get_layer_name();

	// Create the layer in creating mode
	if ( !( LayerManager::CreateAndLockMaskLayer( grid_trans, name, dst_layer,
		this->private_->key_ ) ) )
	{
		dst_layer.reset();
		this->report_error( "Could not allocate enough memory." );
		return false;
	}

	// Record that the layer is locked
	this->private_->created_layers_.push_back( dst_layer );

	// Success
	return true;
}

bool BaseFilter::dispatch_unlock_layer( LayerHandle layer )
{
	bool found_layer = false;
	// Check whether the locked layer is still in the list of layers that this filter locked
	std::vector<LayerHandle>::iterator it;
	
	it = std::find( this->private_->locked_layers_.begin(), 
		this->private_->locked_layers_.end(), layer );

	if ( it != this->private_->locked_layers_.end() )
	{
		// Take the layer out of the list
		this->private_->locked_layers_.erase( it );
		found_layer = true;
	}
	
	// Check whether the locked layer is still in the list of layers that this filter created
	it = std::find( this->private_->created_layers_.begin(), 
		this->private_->created_layers_.end(), layer );

	if ( it != this->private_->created_layers_.end() )
	{
		// Take the layer out of the list
		this->private_->created_layers_.erase( it );
		found_layer = true;
	}

	// If we did not find the layer return false, as we did not lock it in the first place.
	if ( ! found_layer ) 
	{
		this->report_error( "Internal error in filter logic." );
		return false;
	}
	// Send a request to the layer manager to unlock the layer.
	LayerManager::DispatchUnlockLayer( layer, this->private_->key_ );

	// Done
	return true;
}

bool BaseFilter::dispatch_delete_layer( LayerHandle layer )
{
	bool found_layer = false;
	
	// Check whether the locked layer is still in the list of layers that this filter created
	std::vector<LayerHandle>::iterator it = std::find( 
		this->private_->created_layers_.begin(), this->private_->created_layers_.end(), layer );

	if ( it != this->private_->created_layers_.end() )
	{
		// Take the layer out of the list
		this->private_->created_layers_.erase( it );
		found_layer = true;
	}

	it = std::find( this->private_->locked_layers_.begin(), 
		this->private_->locked_layers_.end(), layer );
	if ( it != this->private_->locked_layers_.end() )
	{
		this->private_->locked_layers_.erase( it );
		found_layer = true;
	}
	
	// If we did not find the layer return false, as we did not lock it in the first place.
	if ( ! found_layer ) 
	{
		this->report_error( "Internal error in filter logic." );
		return false;
	}

	// Send a request to the layer manager to unlock the layer.
	LayerManager::DispatchDeleteLayer( layer, this->private_->key_ );

	// Done
	return true;
}


bool BaseFilter::dispatch_insert_data_volume_into_layer( LayerHandle layer, 
	Core::DataVolumeHandle data, bool update_histogram )
{
	// Check whether the layer is of the right type
	DataLayerHandle data_layer = boost::dynamic_pointer_cast<DataLayer>( layer );
	if ( ! data_layer ) return false;

	// Update the data volume if needed.
	// NOTE: We assume that this data volume is not shared by any other thread yet
	// Hence we can update the histogram on the calling thread.
	if ( update_histogram ) data->get_data_block()->update_histogram();
	
	// Ensure that the application thread will process this update.
	LayerManager::DispatchInsertDataVolumeIntoLayer( data_layer, data,
		 this->private_->key_ );
	return true;
}


bool BaseFilter::dispatch_insert_mask_volume_into_layer( LayerHandle layer, 
	Core::MaskVolumeHandle mask )
{	
	// Check whether the layer is of the right type
	MaskLayerHandle mask_layer = boost::dynamic_pointer_cast<MaskLayer>( layer );
	if ( ! mask_layer ) return false;

	// Ensure that the application thread will process this update.
	LayerManager::DispatchInsertMaskVolumeIntoLayer( mask_layer, mask,
		this->private_->key_ );
	return true;
}


void BaseFilter::run()
{
	// Ensure that there not that many processes running in parallel
	BaseFilterLock::Instance()->lock();
	
	try
	{
		this->run_filter();
	}
	catch( ... )
	{
	}
	
	BaseFilterLock::Instance()->unlock();
}

} // end namespace Seg3D

