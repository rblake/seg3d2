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

#include <Core/Geometry/Color.h>
#include <Core/State/StateEngine.h>
#include <Core/State/Actions/ActionAdd.h>
#include <Core/State/Actions/ActionClear.h>
#include <Core/State/Actions/ActionRemove.h>
#include <Core/Viewer/Mouse.h>
#include <Core/Volume/VolumeSlice.h>

#include <Application/LayerManager/LayerManager.h>
#include <Application/Tools/SeedPointsTool.h>
#include <Application/Viewer/Viewer.h>
#include <Application/ViewerManager/ViewerManager.h>

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// Implementation of class 	SeedPointsToolPrivate
//////////////////////////////////////////////////////////////////////////
	
class SeedPointsToolPrivate
{
public:
	Core::VolumeSliceHandle get_target_slice( ViewerHandle viewer );
	void handle_seed_points_changed();
	bool find_point( double world_x, double world_y, 
		Core::VolumeSliceHandle vol_slice, Core::Point& pt );

	SeedPointsTool* tool_;
	ViewerHandle viewer_;
};

Core::VolumeSliceHandle SeedPointsToolPrivate::get_target_slice( ViewerHandle viewer )
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	Core::VolumeSliceHandle vol_slice;
	
	std::string target_layer_id = this->tool_->target_layer_state_->get();
	if ( target_layer_id == Tool::NONE_OPTION_C )
	{
		return vol_slice;
	}
	LayerHandle layer = LayerManager::Instance()->get_layer_by_id( target_layer_id );
	if ( !layer )
	{
		return vol_slice;
	}
	
	switch ( layer->type() )
	{
	case Core::VolumeType::DATA_E:
		vol_slice = viewer->get_data_volume_slice( target_layer_id );
		break;
	case Core::VolumeType::MASK_E:
		vol_slice = viewer_->get_mask_volume_slice( target_layer_id );
		break;
	}
	
	return vol_slice;
}

void SeedPointsToolPrivate::handle_seed_points_changed()
{
	size_t num_of_viewers = ViewerManager::Instance()->number_of_viewers();
	for ( size_t i = 0; i < num_of_viewers; i++ )
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_viewer( i );
		if ( !viewer->is_volume_view() )
		{
			viewer->redraw_overlay();
		}
		else
		{
			viewer->redraw();
		}
	}
}

bool SeedPointsToolPrivate::find_point( double world_x, double world_y, 
									   Core::VolumeSliceHandle vol_slice, Core::Point& pt )
{
	// Step 1. Compute the size of a pixel in world space
	double x0, y0, x1, y1;
	this->viewer_->window_to_world( 0, 0, x0, y0 );
	this->viewer_->window_to_world( 1, 1, x1, y1 );
	double pixel_width = Core::Abs( x1 - x0 );
	double pixel_height = Core::Abs( y1 - y0 );

	// Step 2: Search for the first seed point that's within 4 pixels in each direction
	// from the given search position (world_x, world_y)
	double range_x = pixel_width * 4;
	double range_y = pixel_height * 4;
	std::vector< Core::Point > seed_points;
	{
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		seed_points = this->tool_->seed_points_state_->get();
	}
	size_t num_of_pts = seed_points.size();
	for ( size_t i = 0; i < num_of_pts; i++ )
	{
		double pt_x, pt_y;
		vol_slice->project_onto_slice( seed_points[ i ], pt_x, pt_y );
		if ( Core::Abs( pt_x - world_x ) <= range_x &&
			Core::Abs( pt_y - world_y ) <= range_y )
		{
			pt = seed_points[ i ];
			return true;
		}
	}
	
	return false;
}


//////////////////////////////////////////////////////////////////////////
// Implementation of class SeedPointsTool
//////////////////////////////////////////////////////////////////////////

SeedPointsTool::SeedPointsTool( const std::string& toolid, 
							   size_t version_number, bool auto_number ) :
	Tool( toolid, version_number, auto_number ),
	private_( new SeedPointsToolPrivate )
{
	this->private_->tool_ = this;

	std::vector< LayerIDNamePair > empty_names( 1, 
		std::make_pair( Tool::NONE_OPTION_C, Tool::NONE_OPTION_C ) );
	this->add_state( "target", this->target_layer_state_, Tool::NONE_OPTION_C, empty_names );
	this->add_state( "seed_points", this->seed_points_state_ );

	this->add_connection( this->seed_points_state_->state_changed_signal_.connect(
		boost::bind( &SeedPointsToolPrivate::handle_seed_points_changed, this->private_.get() ) ) );
}

SeedPointsTool::~SeedPointsTool()
{
	this->disconnect_all();
}

bool SeedPointsTool::handle_mouse_enter( size_t viewer_id, int x, int y )
{
	this->private_->viewer_ = ViewerManager::Instance()->get_viewer( viewer_id );
	return true;
}

bool SeedPointsTool::handle_mouse_leave( size_t viewer_id )
{
	this->private_->viewer_.reset();
	return true;
}

bool SeedPointsTool::handle_mouse_press( const Core::MouseHistory& mouse_history, 
										int button, int buttons, int modifiers )
{
	Core::VolumeSliceHandle target_slice;
	double world_x, world_y;

	{
		Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );

		if ( !this->private_->viewer_ || this->private_->viewer_->is_volume_view() )
		{
			return false;
		}
		
		target_slice = this->private_->get_target_slice( this->private_->viewer_ );
		if ( !target_slice || target_slice->out_of_boundary() )
		{
			return false;
		}

		this->private_->viewer_->window_to_world( mouse_history.current_.x_, 
			mouse_history.current_.y_, world_x, world_y );
	}

	if ( modifiers == Core::KeyModifier::NO_MODIFIER_E &&
		button == Core::MouseButton::LEFT_BUTTON_E )
	{
		int u, v;
		target_slice->world_to_index( world_x, world_y, u, v );
		if ( u >= 0 && u < static_cast< int >( target_slice->nx() ) &&
			v >= 0 && v < static_cast< int >( target_slice->ny() ) )
		{
			Core::Point pos;
			target_slice->get_world_coord( world_x, world_y, pos );
			Core::ActionAdd::Dispatch( Core::Interface::GetMouseActionContext(),
				this->seed_points_state_, pos );

			return true;
		}
	}
	else if ( modifiers == Core::KeyModifier::NO_MODIFIER_E &&
		button == Core::MouseButton::RIGHT_BUTTON_E )
	{
		Core::Point pt;
		if ( this->private_->find_point( world_x, world_y, target_slice, pt ) )
		{
			Core::ActionRemove::Dispatch( Core::Interface::GetMouseActionContext(),
				this->seed_points_state_, pt );
		}
		return true;
	}
	
	return false;
}

void SeedPointsTool::redraw( size_t viewer_id, const Core::Matrix& proj_mat )
{
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	if ( viewer->is_volume_view() )
	{
		return;
	}
	Core::VolumeSliceHandle vol_slice = this->private_->get_target_slice( viewer );
	if ( !vol_slice )
	{
		return;
	}
	
	std::vector< Core::Point > seed_points;
	{
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		seed_points = this->seed_points_state_->get();
	}
	size_t num_of_pts = seed_points.size();
	if ( num_of_pts == 0 )
	{
		return;
	}
	
	glPushAttrib( GL_LINE_BIT );
	glLineWidth( 1.0f );

	for ( size_t i = 0; i < num_of_pts; i++ )
	{
		double x_pos, y_pos;
		vol_slice->project_onto_slice( seed_points[ i ], x_pos, y_pos );
		Core::Point pt( x_pos, y_pos, 0 );
		pt = proj_mat * pt;
		int x = static_cast< int >( ( pt[ 0 ] + 1.0 ) * 0.5 * ( viewer->get_width() - 1 ) );
		int y = static_cast< int >( ( pt[ 1 ] + 1.0 ) * 0.5 * ( viewer->get_height() - 1 ) );
		int slice_num = vol_slice->get_closest_slice( seed_points[ i ] );
		bool in_slice = ( !vol_slice->out_of_boundary() && 
			slice_num == static_cast< int >( vol_slice->get_slice_number() ) );
		Core::Color color = in_slice ? Core::Color( 1.0, 1.0, 0.0 ) : Core::Color( 0.6, 0.6, 0.0 );
		glColor3d( color.r(), color.g(), color.b() );
		glBegin( GL_LINES );
			glVertex2i( x - 5, y );
			glVertex2i( x + 5, y );
			glVertex2i( x, y - 5 );
			glVertex2i( x, y + 5 );
		glEnd();
	}
	
	glPopAttrib();
}

void SeedPointsTool::clear( Core::ActionContextHandle context )
{
	Core::ActionClear::Dispatch( context, this->seed_points_state_ );
}

bool SeedPointsTool::has_2d_visual()
{
	return true;
}

bool SeedPointsTool::has_3d_visual()
{
	return true;
}

} // end namespace Seg3D