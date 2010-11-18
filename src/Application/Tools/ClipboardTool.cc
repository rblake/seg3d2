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

#include <Core/Volume/VolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

// Application includes
#include <Application/PreferencesManager/PreferencesManager.h>
#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/ClipboardTool.h>
#include <Application/Tools/Actions/ActionCopy.h>
#include <Application/Tools/Actions/ActionPaste.h>
#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/ViewerManager/ViewerManager.h>
#include <Application/StatusBar/StatusBar.h>

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ClipboardTool )

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// Class ClipboardToolPrivate
//////////////////////////////////////////////////////////////////////////

class ClipboardToolPrivate
{
public:
	void update_slice_numbers();
	void update_slice_type_labels();

	void handle_active_viewer_changed( int active_viewer );
	void handle_viewer_mode_changed( size_t viewer_id, std::string mode );
	void handle_viewer_slice_changed( size_t viewer_id, int slice_num );
	void handle_use_active_viewer_changed( bool use_active_viewer );

	ClipboardTool* tool_;
};

void ClipboardToolPrivate::update_slice_numbers()
{
	ASSERT_IS_APPLICATION_THREAD();

	if ( this->tool_->target_layer_state_->get() == Tool::NONE_OPTION_C )
	{
		return;
	}
	
	Core::VolumeSliceType slice_type( Core::VolumeSliceType::AXIAL_E );
	std::string slice_type_str = this->tool_->slice_type_state_->get();
	if ( slice_type_str == ClipboardTool::CORONAL_C )
	{
		slice_type = Core::VolumeSliceType::CORONAL_E;
	}
	else if ( slice_type_str == ClipboardTool::SAGITTAL_C )
	{
		slice_type = Core::VolumeSliceType::SAGITTAL_E;
	}

	MaskLayerHandle layer = boost::dynamic_pointer_cast< MaskLayer >( 
		LayerManager::Instance()->get_layer_by_id( this->tool_->target_layer_state_->get() ) );
	Core::MaskVolumeSliceHandle vol_slice( new Core::MaskVolumeSlice( 
		layer->get_mask_volume(), slice_type ) );
	this->tool_->copy_slice_number_state_->set_range( 1, static_cast< int >(
		vol_slice->number_of_slices() ) );
	this->tool_->paste_min_slice_number_state_->set_range( 1, static_cast< int >(
		vol_slice->number_of_slices() ) );
	this->tool_->paste_max_slice_number_state_->set_range( 1, static_cast< int >(
		vol_slice->number_of_slices() ) );

	if ( this->tool_->use_active_viewer_state_->get() )
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
		if ( !viewer->is_volume_view() )
		{
			this->tool_->copy_slice_number_state_->set( viewer->slice_number_state_->get() + 1 );
		}
	}
}

void ClipboardToolPrivate::update_slice_type_labels()
{
	Core::OptionLabelPairVector label_options;
	label_options.push_back( std::make_pair( ClipboardTool::AXIAL_C, 
		PreferencesManager::Instance()->z_axis_label_state_->get() ) );
	label_options.push_back( std::make_pair( ClipboardTool::CORONAL_C, 
		PreferencesManager::Instance()->y_axis_label_state_->get() ) );
	label_options.push_back( std::make_pair( ClipboardTool::SAGITTAL_C, 
		PreferencesManager::Instance()->x_axis_label_state_->get() ) );

	this->tool_->slice_type_state_->set_option_list( label_options );
}

void ClipboardToolPrivate::handle_active_viewer_changed( int active_viewer )
{
	if ( !this->tool_->use_active_viewer_state_->get() )
	{
		return;
	}
	
	size_t viewer_id = static_cast< size_t >( active_viewer );
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	if ( !viewer->is_volume_view() )
	{
		this->handle_viewer_mode_changed( viewer_id, viewer->view_mode_state_->get() );
		this->handle_viewer_slice_changed( viewer_id, viewer->slice_number_state_->get() );
	}
}

void ClipboardToolPrivate::handle_viewer_mode_changed( size_t viewer_id, std::string mode )
{
	if ( !this->tool_->use_active_viewer_state_->get() )
	{
		return;
	}
	
	size_t active_viewer = static_cast< size_t >( ViewerManager::Instance()->
		active_viewer_state_->get() );
	if ( viewer_id != active_viewer )
	{
		return;
	}
	
	if ( mode == Viewer::AXIAL_C )
	{
		this->tool_->slice_type_state_->set( ClipboardTool::AXIAL_C );
	}
	else if ( mode == Viewer::CORONAL_C )
	{
		this->tool_->slice_type_state_->set( ClipboardTool::CORONAL_C );
	}
	else if ( mode == Viewer::SAGITTAL_C )
	{
		this->tool_->slice_type_state_->set( ClipboardTool::SAGITTAL_C );
	}
}

void ClipboardToolPrivate::handle_viewer_slice_changed( size_t viewer_id, int slice_num )
{
	if ( !this->tool_->use_active_viewer_state_->get() )
	{
		return;
	}

	size_t active_viewer = static_cast< size_t >( ViewerManager::Instance()->
		active_viewer_state_->get() );
	if ( viewer_id != active_viewer )
	{
		return;
	}

	this->tool_->copy_slice_number_state_->set( slice_num + 1 );
}

void ClipboardToolPrivate::handle_use_active_viewer_changed( bool use_active_viewer )
{
	if ( !use_active_viewer )
	{
		return;
	}
	
	this->tool_->use_active_layer_state_->set( true );
	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	if ( !viewer->is_volume_view() )
	{
		this->handle_viewer_mode_changed( viewer->get_viewer_id(), viewer->view_mode_state_->get() );
		this->handle_viewer_slice_changed( viewer->get_viewer_id(), viewer->slice_number_state_->get() );
	}
}

//////////////////////////////////////////////////////////////////////////
// Class ClipboardTool
//////////////////////////////////////////////////////////////////////////

const std::string ClipboardTool::AXIAL_C( "axial" );
const std::string ClipboardTool::CORONAL_C( "coronal" );
const std::string ClipboardTool::SAGITTAL_C( "sagittal" );


ClipboardTool::ClipboardTool( const std::string& toolid ) :
	SingleTargetTool( Core::VolumeType::MASK_E, toolid ),
	private_( new ClipboardToolPrivate )
{
	this->private_->tool_ = this;
	
	std::string sagittal = SAGITTAL_C + "=" + PreferencesManager::Instance()->x_axis_label_state_->get();
	std::string coronal = CORONAL_C + "=" + PreferencesManager::Instance()->y_axis_label_state_->get();
	std::string axial = AXIAL_C + "=" + PreferencesManager::Instance()->z_axis_label_state_->get();

	this->add_state( "slice_type", this->slice_type_state_, AXIAL_C,  axial + "|" + coronal 
		+ "|" + sagittal );
		
	this->add_connection( PreferencesManager::Instance()->x_axis_label_state_->state_changed_signal_.
		connect( boost::bind( &ClipboardToolPrivate::update_slice_type_labels, this->private_ ) ) );
	this->add_connection( PreferencesManager::Instance()->y_axis_label_state_->state_changed_signal_.
		connect( boost::bind( &ClipboardToolPrivate::update_slice_type_labels, this->private_ ) ) );
	this->add_connection( PreferencesManager::Instance()->z_axis_label_state_->state_changed_signal_.
		connect( boost::bind( &ClipboardToolPrivate::update_slice_type_labels, this->private_ ) ) );	
		
	this->add_state( "copy_slice", this->copy_slice_number_state_, 1, 1, 1, 1 );
	this->add_state( "paste_min_slice", this->paste_min_slice_number_state_, 1, 1, 1, 1 );
	this->add_state( "paste_max_slice", this->paste_max_slice_number_state_, 1, 1, 1, 1 );
	this->add_state( "use_active_viewer", this->use_active_viewer_state_, true );

	this->private_->update_slice_numbers();
	this->private_->handle_use_active_viewer_changed( true );

	this->add_connection( this->target_layer_state_->state_changed_signal_.connect(
		boost::bind( &ClipboardToolPrivate::update_slice_numbers, this->private_ ) ) );
	this->add_connection( this->slice_type_state_->state_changed_signal_.connect(
		boost::bind( &ClipboardToolPrivate::update_slice_numbers, this->private_ ) ) );
	this->add_connection( this->use_active_viewer_state_->value_changed_signal_.connect( 
		boost::bind( &ClipboardToolPrivate::handle_use_active_viewer_changed, this->private_, _1 ) ) );
	this->add_connection( ViewerManager::Instance()->active_viewer_state_->value_changed_signal_.
		connect( boost::bind( &ClipboardToolPrivate::handle_active_viewer_changed, 
		this->private_, _1 ) ) );
	size_t num_of_viewrs = ViewerManager::Instance()->number_of_viewers();
	for ( size_t i = 0; i < num_of_viewrs; ++i )
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_viewer( i );
		this->add_connection( viewer->view_mode_state_->value_changed_signal_.connect( 
			boost::bind( &ClipboardToolPrivate::handle_viewer_mode_changed, 
			this->private_, i, _2 ) ) );
		this->add_connection( viewer->slice_number_state_->value_changed_signal_.connect(
			boost::bind( &ClipboardToolPrivate::handle_viewer_slice_changed,
			this->private_, i, _1 ) ) );
	}
}

ClipboardTool::~ClipboardTool()
{
	this->disconnect_all();
}

void ClipboardTool::copy( Core::ActionContextHandle context )
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	const std::string& layer_id = this->target_layer_state_->get();
	if ( layer_id == Tool::NONE_OPTION_C )
	{
		return;
	}
	
	ActionCopy::Dispatch( context, layer_id, this->slice_type_state_->index(),
		static_cast< size_t >( this->copy_slice_number_state_->get() - 1 ) );
}

void ClipboardTool::paste( Core::ActionContextHandle context )
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	const std::string& layer_id = this->target_layer_state_->get();
	if ( layer_id == Tool::NONE_OPTION_C )
	{
		return;
	}

	ActionPaste::Dispatch( context, layer_id, this->slice_type_state_->index(),
		static_cast< size_t >( this->paste_min_slice_number_state_->get() - 1 ),
		static_cast< size_t >( this->paste_max_slice_number_state_->get() - 1 ) );
}

void ClipboardTool::grab_min_paste_slice()
{
	if ( !Core::Application::IsApplicationThread() )
	{
		Core::Application::PostEvent( boost::bind( &ClipboardTool::grab_min_paste_slice, this ) );
		return;
	}
	
	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	if ( viewer->is_volume_view() )
	{
		StatusBar::SetMessage( Core::LogMessageType::ERROR_E, "Active viewer not in 2D view." );
		return;
	}
	
	this->paste_min_slice_number_state_->set( viewer->slice_number_state_->get() + 1 );
}

void ClipboardTool::grab_max_paste_slice()
{
	if ( !Core::Application::IsApplicationThread() )
	{
		Core::Application::PostEvent( boost::bind( &ClipboardTool::grab_max_paste_slice, this ) );
		return;
	}

	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	if ( viewer->is_volume_view() )
	{
		StatusBar::SetMessage( Core::LogMessageType::ERROR_E, "Active viewer not in 2D view." );
		return;
	}

	this->paste_max_slice_number_state_->set( viewer->slice_number_state_->get() + 1 );
}

} // end namespace Seg3D
