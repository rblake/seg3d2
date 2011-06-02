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

// Application includes
#include <Application/Tool/ToolFactory.h>
#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/Layer/LayerGroup.h>

// StateEngine of the tool
#include <Application/Tools/ICPFilter.h>

// Action associated with tool
#include <Application/Filters/Actions/ActionICPRegisterFilter.h>
#include <Application/Filters/Actions/ActionICPTransformFilter.h>

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ICPFilter )

namespace Seg3D
{

ICPFilter::ICPFilter( const std::string& toolid ) :
	SingleTargetTool( Core::VolumeType::MASK_E, toolid )
{
	// Create an empty list of label options
	std::vector< LayerIDNamePair > empty_list( 1, 
		std::make_pair( Tool::NONE_OPTION_C, Tool::NONE_OPTION_C ) );

	// Need to set ranges and default values for all parameters
	//this->add_state( "replace", this->replace_state_, false );

	// Whether we use a mask to find which components to use
	this->add_state( "mask", this->mask_state_, Tool::NONE_OPTION_C, empty_list );
	//this->add_dependent_layer_input( this->mask_state_, Core::VolumeType::MASK_E, true );
	this->add_dependent_layer_input( this->mask_state_, Core::VolumeType::MASK_E, true, true );

	this->add_state( "iterations", this->iterations_state_, 5, 1, 1000, 1 );

	std::vector< std::string > empty_option;
	this->add_state( "target_layers", this->target_layers_state_, empty_option, "" );

	this->add_state( "registration_ready", this->registration_ready_state_, false );

	this->add_connection( this->mask_state_->value_changed_signal_.connect(
		boost::bind( &ICPFilter::handle_mask_layer_changed, this, _2 ) ) );

	this->add_connection( this->target_layer_state_->value_changed_signal_.connect(
		boost::bind( &ICPFilter::handle_target_layer_changed, this, _2 ) ) );

	this->add_connection( this->iterations_state_->value_changed_signal_.connect(
		boost::bind( &ICPFilter::handle_iteration_changed, this ) ) );

	this->add_state( "transformation_matrix", this->transform_matrix_, std::vector<double>());

	this->add_connection( LayerManager::Instance()->layers_changed_signal_.connect(
		boost::bind( &ICPFilter::handle_layers_changed, this ) ) );

}
	
ICPFilter::~ICPFilter()
{
	this->disconnect_all();
}

void ICPFilter::handle_target_layer_changed( std::string layer_id )
{
	registration_ready_state_->set( false );
}


void ICPFilter::handle_iteration_changed(  )
{
	registration_ready_state_->set( false );
}

void ICPFilter::handle_layers_changed()
{
	std::string mask_layer_id = this->mask_state_->get();
	
	LayerHandle layer = LayerManager::Instance()->get_layer_by_id( mask_layer_id );

	if ( layer )
	{
		LayerGroupHandle group = layer->get_layer_group();
		std::string group_id = group->get_group_id();

		Core::VolumeType mask_type = layer->get_type();

		std::vector< LayerIDNamePair > layer_names;
		std::vector< std::string > selected_layers;
		if ( group_id != "" && group_id != Tool::NONE_OPTION_C )
		{
			group->get_layer_names( layer_names, mask_type );

			LayerHandle target_layer = LayerManager::Instance()->get_layer_by_id( 
				this->target_layer_state_->get());

			std::string target_layer_name = target_layer->get_layer_name();
			std::vector<LayerIDNamePair>::iterator it = layer_names.begin();

			for ( ; it!=layer_names.end(); ++it )
			{
				if ( target_layer_name.compare( (*it).second ) == 0 ) 
				{
					break;
				}
				//std::find( layer_names.begin(), layer_names.end(), target_layer_name );
			}

			if ( it!= layer_names.end() )
			{
				layer_names.erase( it );
			}
		}

		this->target_layers_state_->set_option_list( layer_names );
		if ( selected_layers.size() > 0 )
		{
			this->target_layers_state_->set( selected_layers );
		}
	}

}

void ICPFilter::handle_mask_layer_changed( std::string layer_id )
{
	LayerHandle layer = LayerManager::Instance()->get_layer_by_id( layer_id );

	if ( layer )
	{
		registration_ready_state_->set( false );

		LayerGroupHandle group = layer->get_layer_group();
		std::string group_id = group->get_group_id();

		Core::VolumeType mask_type = layer->get_type();

		std::vector< LayerIDNamePair > layer_names;
		std::vector< std::string > selected_layers;
		if ( group_id != "" && group_id != Tool::NONE_OPTION_C )
		{
			group->get_layer_names( layer_names, mask_type );

			LayerHandle target_layer = LayerManager::Instance()->get_layer_by_id( 
				this->target_layer_state_->get());

			std::string target_layer_name = target_layer->get_layer_name();
			std::vector<LayerIDNamePair>::iterator it = layer_names.begin();

			for ( ; it!=layer_names.end(); ++it )
			{
				if ( target_layer_name.compare( (*it).second ) == 0 ) 
				{
					break;
				}
			}

			if ( it!= layer_names.end() )
			{
				layer_names.erase( it );
			}
		}

		this->target_layers_state_->set_option_list( layer_names );
		if ( selected_layers.size() > 0 )
		{
			this->target_layers_state_->set( selected_layers );
		}
	}

}

void ICPFilter::execute( Core::ActionContextHandle context )
{
	// NOTE: Need to lock state engine as this function is run from the interface thread
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	ActionICPRegisterFilter::Dispatch( context,
		this->target_layer_state_->get(),
		this->mask_state_->get(),
		this->iterations_state_->get(),
		this->toolid() );	

}

void ICPFilter::apply( Core::ActionContextHandle context )
{
	// NOTE: Need to lock state engine as this function is run from the interface thread
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	ActionICPTransformFilter::Dispatch( context,
		this->target_layer_state_->get(),
		this->target_layers_state_->get(),
		this->toolid()
		);		
}

} // end namespace Seg3D


