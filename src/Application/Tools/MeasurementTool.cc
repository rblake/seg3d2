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

// Boost includes
#include <boost/foreach.hpp>
#include <boost/regex.hpp>

// Application includes
#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/Actions/ActionSetMeasurementPoint.h>
#include <Application/Tools/Actions/ActionSetMeasurementVisible.h>
#include <Application/Tools/MeasurementTool.h>
#include <Application/ViewerManager/ViewerManager.h>

// Core includes
#include <Core/Geometry/Point.h>
#include <Core/RenderResources/RenderResources.h>
#include <Core/State/Actions/ActionAdd.h>
#include <Core/State/Actions/ActionRemove.h>
#include <Core/State/Actions/ActionSetAt.h>
#include <Core/TextRenderer/TextRenderer.h>
#include <Core/Utils/Lockable.h>
#include <Core/Viewer/Mouse.h>

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, MeasurementTool )

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// Class MeasurementPoint
//////////////////////////////////////////////////////////////////////////

class MeasurementPoint
{
public:
	MeasurementPoint();

	bool is_valid();
	void invalidate();

	// Storing id instead of index just in case measurements change on another thread
	std::string measurement_id_;
	int point_index_;
};

MeasurementPoint::MeasurementPoint() 
{
	this->invalidate();
}

bool MeasurementPoint::is_valid()
{
	return ( this->measurement_id_ != "" && ( this->point_index_ == 0 || this->point_index_ == 1 ) );	
}

void MeasurementPoint::invalidate()
{
	this->measurement_id_ = "";
	this->point_index_ = -1;
}

//////////////////////////////////////////////////////////////////////////
// Class MeasurementToolPrivate
//////////////////////////////////////////////////////////////////////////

class MeasurementToolPrivate : public Core::RecursiveLockable
{
public:
	//
	// Tool Interface 
	// 

	void handle_measurements_changed();
	void handle_units_selection_changed( std::string units );
	void handle_active_layer_changed( LayerHandle active_layer );
	void handle_opacity_changed();
	void handle_active_viewer_changed( int active_viewer );
	void handle_viewer_slice_changed( size_t viewer_id, int slice_num );

	void initialize_id_counter();
	std::string get_next_measurement_id();

	//
	// Viewer 
	//

	// UPDATE_VIEWERS:
	// Redraw in relevant viewers.
	void update_viewers();

	// FIND_HOVER_POINT:
	// Find measurement point (if any) that mouse is currently hovering over.
	bool find_hover_point();

	// MOVE_HOVER_POINT_TO_MOUSE:
	// Move hover point to current mouse position.
	void move_hover_point_to_mouse();

	// SNAP_HOVER_POINT_TO_SLICE:
	// Snap the hover point to its projected position on the current slice.
	bool snap_hover_point_to_slice();

	// UPDATE_HOVER_POINT:
	// Find or move the hover point, depending on whether user is editing or not.
	void update_hover_point();

	// Locks: StateEngine
	bool get_hover_measurement( size_t& index, Core::Measurement& measurement, 
		Core::Point& world_point ); 

	void start_editing();
	void finish_editing();

	bool get_mouse_world_point( Core::Point& world_point );

	// IN_SLICE:
	// Return true if the point is in the current slice.  This is slightly different from seed 
	// points -- in this case the point needs to be within epsilon of the slice, not just nearest
	// to it.
	bool in_slice( ViewerHandle viewer, const Core::Point& world_point );

	MeasurementTool* tool_;

	ViewerHandle viewer_; // Should be mutex-protected

	std::string active_group_id_; // Only accessed from application thread

	// Is a measurement being edited?
	bool editing_; // Should be mutex-protected

	// Measurement point that is currently hovered over by the mouse. 
	// We need to know the hover point for left-click + mouse move (editing),
	// middle-click (move point to current slice) and right-click (delete measurement).
	// Want to store this so that we don't have to "find" the hover point every time the mouse moves
	// since during editing we might find a different point and start moving that one instead.  
	// Also need to check this in redraw when determining which measurement (if any) to draw as a 
	// dotted line.
	MeasurementPoint hover_point_; // Should be mutex-protected

	int measurement_id_counter_; // Should be mutex-protected

	Core::MousePosition mouse_pos_; // Should be mutex-protected

	int saved_num_measurements_; // Only accessed from application thread

	bool handle_measurements_changed_blocked_; // Only accessed from application thread
};

// Called in response to state changed signal
void MeasurementToolPrivate::handle_measurements_changed()
{
	// Running on application thread, so measurements list and active index won't be changed out 
	// from under us.
	ASSERT_IS_APPLICATION_THREAD();

	// Prevent circular updates
	if( this->handle_measurements_changed_blocked_ ) 
	{
		return;
	}
	this->handle_measurements_changed_blocked_ = true;

	int num_measurements = static_cast< int >( this->tool_->measurements_state_->get().size() );
	
	bool num_measurements_changed = num_measurements != saved_num_measurements_;
	if( num_measurements_changed )
	{
		this->tool_->num_measurements_changed_signal_();
	}

	// If measurements were deleted, renumber IDs
	if( num_measurements < saved_num_measurements_ )
	{
		lock_type lock( this->get_mutex() );
		this->measurement_id_counter_ = 0;

		if( num_measurements > 0 )
		{
			std::vector< Core::Measurement > measurements = 
				this->tool_->measurements_state_->get();

			for( size_t m_idx = 0; m_idx < measurements.size(); m_idx++ )
			{
				std::string id = "M" + Core::ExportToString( this->measurement_id_counter_ );
				measurements[ m_idx ].set_id( id );
				this->measurement_id_counter_++;
			}
			// Do all measurement changes at once to minimize updates
			this->tool_->measurements_state_->set( measurements );	
		}
	}

	// Measurements may have been added or removed, so update the active index
	if( num_measurements > 0 )
	{
		int active_index = this->tool_->active_index_state_->get();
		bool active_index_invalid = active_index == -1 || active_index >= num_measurements;

		if( num_measurements_changed || active_index_invalid )
		{
			// Set active index to end of list.	
			Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(), 
				this->tool_->active_index_state_, num_measurements - 1 );
		}
	}
	else
	{
		Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(), 
			this->tool_->active_index_state_, -1 );
	}

	this->saved_num_measurements_ = num_measurements;
	this->handle_measurements_changed_blocked_ = false;

	// Need to redraw the overlay
	this->update_viewers();
}

void MeasurementToolPrivate::handle_units_selection_changed( std::string units )
{
	// Don't need to lock state engine because we're running on app thread
	ASSERT_IS_APPLICATION_THREAD();

	bool old_show_world_units_state = this->tool_->show_world_units_state_->get();

	if ( units == MeasurementTool::WORLD_UNITS_C )
	{
		this->tool_->show_world_units_state_->set( true );
	}
	else
	{
		this->tool_->show_world_units_state_->set( false );
	}

	// If units have changed, emit signal
	if( old_show_world_units_state != this->tool_->show_world_units_state_->get() )
	{
		this->tool_->units_changed_signal_();

		// Need to redraw the overlay since length is rendered with measurements
		this->update_viewers();
	}
}

void MeasurementToolPrivate::handle_active_layer_changed( LayerHandle active_layer )
{
	ASSERT_IS_APPLICATION_THREAD();

	// To minimize measurement table updates, only emit units_changed_signal_ when 
	// the active group has changed AND the index units are selected.
	bool show_index_units = !this->tool_->show_world_units_state_->get();
	if ( active_layer )
	{
		// If the active group has changed
		std::string curr_active_group_id = active_layer->get_layer_group()->get_group_id();
		// Don't need to mutex-protect active_group_id_ because it will only be accessed
		// on the application thread.
		// Check to see if the group size is 1 to handle case where layer is deleted, then 
		// delete is undone but group id didn't get changed in process.
		if( this->active_group_id_ != curr_active_group_id || 
			active_layer->get_layer_group()->get_list_size() == 1 )  
		{
			this->active_group_id_ = curr_active_group_id;

			if( show_index_units )
			{
				// Signal that measurement units have changed
				this->tool_->units_changed_signal_();
			}	
		}
	}
	else
	{
		this->active_group_id_ = "";
	}
}

void MeasurementToolPrivate::handle_opacity_changed()
{
	ASSERT_IS_APPLICATION_THREAD();

	this->update_viewers();
}

void MeasurementToolPrivate::handle_active_viewer_changed( int active_viewer )
{
	ASSERT_IS_APPLICATION_THREAD();

	size_t viewer_id = static_cast< size_t >( active_viewer );
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	if ( viewer && !viewer->is_volume_view() )
	{
		this->handle_viewer_slice_changed( viewer_id, viewer->slice_number_state_->get() );
	}
}

void MeasurementToolPrivate::handle_viewer_slice_changed( size_t viewer_id, int slice_num )
{
	ASSERT_IS_APPLICATION_THREAD();

	size_t active_viewer = static_cast< size_t >( ViewerManager::Instance()->
		active_viewer_state_->get() );
	if ( viewer_id != active_viewer )
	{
		return;
	}

	// Update hover point based on new slice
	this->update_hover_point();
}

void MeasurementToolPrivate::initialize_id_counter()
{
	ASSERT_IS_INTERFACE_THREAD();

	lock_type lock( this->get_mutex() );

	// Find highest ID in use
	this->measurement_id_counter_ = 0;
	Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
	const std::vector< Core::Measurement >& measurements = this->tool_->measurements_state_->get();
	BOOST_FOREACH( Core::Measurement m, measurements )
	{
		int measurement_id = 0;
		boost::regex reg( "(M)(\\d*)" );
		boost::smatch sub_matches;
		// Note that we can't just pass m.get_id() directly into regex_match because the string will
		// go out of scope and sub_matches contains iterators that point into the string.  
		std::string str = m.get_id();
		if( boost::regex_match( str, sub_matches, reg ) ) 
		{
			Core::ImportFromString( sub_matches[ 2 ].str(), measurement_id );

			if( measurement_id >= this->measurement_id_counter_ )
			{
				this->measurement_id_counter_ = measurement_id + 1;
			}
		}
	}
}

std::string MeasurementToolPrivate::get_next_measurement_id()
{
	ASSERT_IS_INTERFACE_THREAD();

	lock_type lock( this->get_mutex() );

	if( this->measurement_id_counter_ == -1 )
	{
		this->initialize_id_counter();
	}

	// Find first id not already in use (may have been loaded from session file)
	while( true )
	{
		std::string canditate_id = "M" + 
			Core::ExportToString( this->measurement_id_counter_ );
		this->measurement_id_counter_++;

		// Make sure this ID isn't in use
		bool in_use = false;
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		const std::vector< Core::Measurement >& measurements = 
			this->tool_->measurements_state_->get();
		BOOST_FOREACH( Core::Measurement m, measurements )
		{
			if( canditate_id == m.get_id() )
			{
				in_use = true;
				break;
			}
		}
		if( !in_use )
		{
			return canditate_id;
		}
	}
}

void MeasurementToolPrivate::update_viewers()
{
	// May be called from application or interface thread
	ViewerManager::Instance()->update_2d_viewers_overlay();
}

bool MeasurementToolPrivate::find_hover_point()
{
	// May be called from application (slice changed) or interface (mouse move) thread
	lock_type lock( this->get_mutex() );

	// Need to find first point within radius.
	// NOTE: Should we find closest point within radius instead?

	if( !this->viewer_ ) 
	{
		this->hover_point_.invalidate();
		return false;
	}

	Core::VolumeSliceHandle volume_slice = this->viewer_->get_active_volume_slice();
	if( !volume_slice )
	{
		this->hover_point_.invalidate();
		return false;
	}

	// Compute the size of a pixel in world space
	double x0, y0, x1, y1;
	this->viewer_->window_to_world( 0, 0, x0, y0 );
	this->viewer_->window_to_world( 1, 1, x1, y1 );
	double pixel_width = Core::Abs( x1 - x0 );
	double pixel_height = Core::Abs( y1 - y0 );

	// Compute the mouse position in world space
	double world_x, world_y;
	this->viewer_->window_to_world( this->mouse_pos_.x_, this->mouse_pos_.y_, world_x, world_y );

	// Search for the first vertex that's within 2 pixels of current mouse position
	double range_x = pixel_width * 4;
	double range_y = pixel_height * 4;

	Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
	const std::vector< Core::Measurement >& measurements = this->tool_->measurements_state_->get();
	BOOST_FOREACH( Core::Measurement m, measurements )
	{
		// Ignore hidden measurements
		if( m.get_visible() )
		{
			// For each measurement point, project onto slice and check to see if in radius
			for( int point_index = 0; point_index < 2; point_index++ )
			{
				Core::Point measurement_point;
				m.get_point( point_index, measurement_point );
				double pt_x, pt_y;
				volume_slice->project_onto_slice( measurement_point, pt_x, pt_y );
				if ( Core::Abs( pt_x - world_x ) <= range_x &&
					Core::Abs( pt_y - world_y ) <= range_y )
				{
					this->hover_point_.measurement_id_ = m.get_id();
					this->hover_point_.point_index_ = point_index; 
					return true;
				}
			}
		}
	}

	// If no point hovered over, invalidate hover point
	this->hover_point_.invalidate();
	return false;
}

void MeasurementToolPrivate::move_hover_point_to_mouse()
{
	// May be called from application (slice changed) or interface (mouse move) thread
	lock_type lock( this->get_mutex() );
	
	if( this->hover_point_.is_valid() )
	{
		// Convert current mouse coords to 3D world point
		Core::Point moved_pt;
		if( this->get_mouse_world_point( moved_pt ) )
		{
			// Update hover point
			ActionSetMeasurementPoint::Dispatch( Core::Interface::GetMouseActionContext(), 
				this->tool_->measurements_state_, this->hover_point_.measurement_id_, 
				this->hover_point_.point_index_, moved_pt );
		}
	}
}

bool MeasurementToolPrivate::snap_hover_point_to_slice()
{
	ASSERT_IS_INTERFACE_THREAD();
	lock_type lock( this->get_mutex() );

	// Find world hover point
	size_t measurement_index;
	Core::Measurement edited_measurement;
	Core::Point measurement_point;
	if( this->get_hover_measurement( measurement_index, edited_measurement, measurement_point ) )
	{
		if( !this->viewer_ )
		{
			return false;
		}

		// Project hover point onto current slice, find world point
		Core::VolumeSliceHandle active_slice = this->viewer_->get_active_volume_slice();
		if( !active_slice )
		{
			return false;
		}
		
		double world_x, world_y;
		active_slice->project_onto_slice( measurement_point, world_x, world_y );
		Core::Point moved_pt;
		active_slice->get_world_coord( world_x, world_y, moved_pt );

		// Update hover point
		ActionSetMeasurementPoint::Dispatch( Core::Interface::GetMouseActionContext(), 
			this->tool_->measurements_state_, this->hover_point_.measurement_id_, 
			this->hover_point_.point_index_, moved_pt );
		return true;
	}
	return false;
}

void MeasurementToolPrivate::update_hover_point()
{
	// May be called from application (slice changed) or interface (mouse move) thread
	lock_type lock( this->get_mutex() );

	if( !this->viewer_ )
	{
		return;
	}

	if( this->editing_ )
	{
		// Move hovered point.  We don't want to "find" the hover point every time the mouse moves
		// during editing since we might find a different point and start moving that one instead.
		this->move_hover_point_to_mouse();
	}
	else
	{	
		// Find the measurement point that the mouse is currently hovering over, if there is one
		if( this->find_hover_point() )
		{
			// Set cursor to open hand to indicate that point could be selected (but isn�t)
			this->viewer_->set_cursor( Core::CursorShape::OPEN_HAND_E );
		}
		else
		{
			// Set cross icon color to grey (normal)
			this->viewer_->set_cursor( Core::CursorShape::CROSS_E );
		}
	}
}

bool MeasurementToolPrivate::get_hover_measurement( size_t& index, Core::Measurement& measurement, 
	Core::Point& world_point )
{
	// May be called from application or interface thread
	lock_type lock( this->get_mutex() );

	if( this->hover_point_.is_valid() )
	{
		// Get measurements
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		const std::vector< Core::Measurement >& measurements = 
			this->tool_->measurements_state_->get();

		// Find one with matching id
		for( size_t i = 0; i < measurements.size(); i++ )
		{
			if( measurements[ i ].get_id() == this->hover_point_.measurement_id_ )
			{
				index = i;
				measurement = measurements[ i ];
				measurement.get_point( this->hover_point_.point_index_, world_point );
				return true;
			}
		}
	}
	return false;
}

void MeasurementToolPrivate::start_editing()
{
	ASSERT_IS_INTERFACE_THREAD();
	lock_type lock( this->get_mutex() );

	if( !this->viewer_ )
	{
		return;
	}

	this->editing_ = true;
	
	// Use blank cursor so that features of interest are not obscured
	this->viewer_->set_cursor( Core::CursorShape::BLANK_E );

	// Redraw so that measurement is drawn with dotted line instead of solid line
	this->update_viewers();
}

void MeasurementToolPrivate::finish_editing()
{
	ASSERT_IS_INTERFACE_THREAD();
	lock_type lock( this->get_mutex() );

	if( !this->viewer_ )
	{
		return;
	}

	this->editing_ = false;

	// No need to update measurement (interactively updated) 
	// Change cursor to indicate that editing is not happening
	// NOTE: Cursor has three states: editing, hovering, and normal	
	if( this->hover_point_.is_valid() )
	{
		this->viewer_->set_cursor( Core::CursorShape::OPEN_HAND_E );
	}
	else
	{
		this->viewer_->set_cursor( Core::CursorShape::CROSS_E );
	}
	
	// Redraw so that measurement is drawn with solid line instead of dotted line
	this->update_viewers();
}

bool MeasurementToolPrivate::get_mouse_world_point( Core::Point& world_point )
{
	// May be called from application or interface thread
	lock_type lock( this->get_mutex() );

	if( !this->viewer_ )
	{
		return false;
	}

	// Make sure to use current viewer, slice, and mouse coordinates.  These aren't passed
	// as parameters because we don't have this info in the case where the view has been changed
	// by a key command.
	double world_x, world_y;
	this->viewer_->window_to_world( this->mouse_pos_.x_, this->mouse_pos_.y_, world_x, world_y );
	Core::VolumeSliceHandle active_slice = this->viewer_->get_active_volume_slice();
	if( !active_slice )
	{
		return false;
	}

	active_slice->get_world_coord( world_x, world_y, world_point );	
	return true;
}

bool MeasurementToolPrivate::in_slice( ViewerHandle viewer, const Core::Point& world_point )
{
	// May be called from interface thread or rendering thread

	// Basically in_slice has to be within epsilon of this slice so that editing won't move the
	// measurement points
	Core::VolumeSliceHandle volume_slice = viewer->get_active_volume_slice();
	if( !volume_slice )
	{
		return false;
	}

	double slice_depth = volume_slice->depth();
	double point_depth;
	double i_pos, j_pos;
	volume_slice->project_onto_slice( world_point, i_pos, j_pos, point_depth );

	if( point_depth == slice_depth ) 
	{
		return true;
	}
	return false;
}

//void create_test_data( std::vector< Core::Measurement >& measurements )
//{
//	// Populate measurements list with test data
//	measurements.push_back( 
//		Core::Measurement( true, "M1", "Knee", Core::Point(0, 0, 0), 
//		Core::Point(1, 1, 1), Core::AXIAL_E, 50, 1 ) );
//	measurements.push_back( 
//		Core::Measurement( true, "M2", "Heart", Core::Point(0, 0, 0), 
//		Core::Point(2, 2, 2), Core::AXIAL_E, 50, 1 ) );
//	measurements.push_back( 
//		Core::Measurement( true, "M3", "Head", Core::Point(0, 0, 0), 
//		Core::Point(3, 3, 3), Core::AXIAL_E, 50, 1 ) );	
//	measurements.push_back( 
//		Core::Measurement( true, "M4", "Toe", Core::Point(0, 0, 0), 
//		Core::Point(4, 4, 4), Core::AXIAL_E, 50, 1 ) );	
//	measurements.push_back( 
//		Core::Measurement( true, "M5", "Eye", Core::Point(0, 0, 0), 
//		Core::Point(5, 5, 5), Core::AXIAL_E, 50, 1 ) );	
//	measurements.push_back( 
//		Core::Measurement( true, "M6", "Nose", Core::Point(0, 0, 0), 
//		Core::Point(6, 6, 6), Core::AXIAL_E, 50, 1 ) );	
//	measurements.push_back( 
//		Core::Measurement( true, "M7", "Hand", Core::Point(0, 0, 0), 
//		Core::Point(7, 7, 7), Core::AXIAL_E, 50, 1 ) );	
//	measurements.push_back( 
//		Core::Measurement( true, "M8", "Ear", Core::Point(0, 0, 0), 
//		Core::Point(8, 8, 8), Core::AXIAL_E, 50, 1 ) );	
//}

//////////////////////////////////////////////////////////////////////////
// Class MeasurementTool
//////////////////////////////////////////////////////////////////////////

const std::string MeasurementTool::INDEX_UNITS_C( "index_units" );
const std::string MeasurementTool::WORLD_UNITS_C( "world_units" );

MeasurementTool::MeasurementTool( const std::string& toolid ) :
	Tool( toolid ),
	private_( new MeasurementToolPrivate )
{
	this->private_->tool_ = this;
	this->private_->active_group_id_ = "";
	this->private_->editing_ = false;
	this->private_->measurement_id_counter_ = -1;
	this->private_->saved_num_measurements_ = 0;
	this->private_->handle_measurements_changed_blocked_ = false;
	
	// State variable gets allocated here
	this->add_state( "measurements", this->measurements_state_ );
	this->add_state( "active_index", this->active_index_state_, -1 );
	this->add_state( "units_selection", this->units_selection_state_, WORLD_UNITS_C, 
		INDEX_UNITS_C + "=Index Units|" +
		WORLD_UNITS_C + "=World Units" );
	this->add_state( "show_world_units", this->show_world_units_state_, true );
	this->add_state( "opacity", this->opacity_state_, 1.0, 0.0, 1.0, 0.1 );

	LayerHandle active_layer = LayerManager::Instance()->get_active_layer();
	if( active_layer )
	{
		this->private_->active_group_id_ = active_layer->get_layer_group()->get_group_id();
	}

	this->add_connection( this->measurements_state_->state_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_measurements_changed, this->private_ ) ) );
	this->add_connection( this->active_index_state_->state_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::update_viewers, this->private_ ) ) );
	this->add_connection( this->units_selection_state_->value_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_units_selection_changed, this->private_, _2 ) ) );
	this->add_connection( LayerManager::Instance()->active_layer_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_active_layer_changed, this->private_, _1 ) ) );
	this->add_connection( this->opacity_state_->state_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_opacity_changed, this->private_ ) ) );
	this->add_connection( ViewerManager::Instance()->active_viewer_state_->value_changed_signal_.
		connect( boost::bind( &MeasurementToolPrivate::handle_active_viewer_changed, 
		this->private_, _1 ) ) );

	size_t num_of_viewers = ViewerManager::Instance()->number_of_viewers();
	for ( size_t i = 0; i < num_of_viewers; ++i )
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_viewer( i );
		this->add_connection( viewer->slice_number_state_->value_changed_signal_.connect(
			boost::bind( &MeasurementToolPrivate::handle_viewer_slice_changed,
			this->private_, i, _1 ) ) );
	}
}

MeasurementTool::~MeasurementTool()
{
	this->disconnect_all();
}

bool MeasurementTool::handle_mouse_move( ViewerHandle viewer, 
	const Core::MouseHistory& mouse_history, int button, int buttons, int modifiers )
{
	if( viewer->is_volume_view() ) 
	{
		return false;
	}

	MeasurementToolPrivate::lock_type lock( this->private_->get_mutex() );

	this->private_->viewer_ = viewer;
	this->private_->mouse_pos_.x_ = mouse_history.current_.x_;
	this->private_->mouse_pos_.y_ = mouse_history.current_.y_;
	this->private_->update_hover_point();

	// May have moved into new viewer, need to set its cursor if editing
	if( this->private_->editing_ && this->private_->viewer_ )
	{
		// Use blank cursor so that features of interest are not obscured
		this->private_->viewer_->set_cursor( Core::CursorShape::BLANK_E );
	}

	// Pass handling on to normal handler 
	return false;
}

bool MeasurementTool::handle_mouse_press( ViewerHandle viewer, 
	const Core::MouseHistory& mouse_history, int button, int buttons, int modifiers )
{
	if( viewer->is_volume_view() ) 
	{
		return false;
	}

	MeasurementToolPrivate::lock_type lock( this->private_->get_mutex() );
	this->private_->viewer_ = viewer;

	if( modifiers != Core::KeyModifier::NO_MODIFIER_E ) 
	{
		return false;
	}

	if( button == Core::MouseButton::LEFT_BUTTON_E )
	{
		if( this->private_->editing_ ) 
		{
			// State: We were editing

			// Done editing
			this->private_->finish_editing();

		}
		else 
		{
			// Find 3D world hover point
			size_t measurement_index;
			Core::Measurement measurement;
			Core::Point hover_point;
			if( this->private_->get_hover_measurement( measurement_index, measurement, hover_point ) 
				&& this->private_->in_slice( viewer, hover_point ) )
			{
				// State: Hovering over point that is in slice 
				
				// Make hover measurement active
				Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(), 
					this->active_index_state_, static_cast< int >( measurement_index ) );

				// Start editing
				this->private_->start_editing();
			}
			else 
			{
				// State: Hovering over no point or point not in slice

				// Get 3D point from 2D mouse coords, create new measurement with P1 = P2, add to state vector using action.
				Core::Point pt;
				if( this->private_->get_mouse_world_point( pt ) )
				{
					// Create new measurement 
					Core::Measurement measurement;

					// Need ID for measurement
					measurement.set_id( this->private_->get_next_measurement_id() );
					measurement.set_point( 0, pt );
					measurement.set_point( 1, pt );
					measurement.set_visible( true );

					// Add measurement to state vector
					Core::ActionAdd::Dispatch( Core::Interface::GetMouseActionContext(), 
						this->measurements_state_, measurement );

					// Second point in measurement needs to be the hover point
					this->private_->hover_point_.measurement_id_ = measurement.get_id();
					this->private_->hover_point_.point_index_ = 1; // Editing second point

					this->private_->start_editing();
				}
			}
		}
		return true;
	}
	else if( button == Core::MouseButton::MID_BUTTON_E )
	{
		if( this->private_->hover_point_.is_valid() ) // Hovering
		{
			// This point does not have to be in current slice
			// Snap point to current slice
			this->private_->snap_hover_point_to_slice();

			// TODO: Different icon for hovering but not in slice?  
		}
	}
	else if( button == Core::MouseButton::RIGHT_BUTTON_E )
	{
		// Find 3D world hover point
		size_t measurement_index;
		Core::Measurement measurement;
		Core::Point hover_point;
		if( this->private_->get_hover_measurement( measurement_index, measurement, hover_point ) && 
			this->private_->in_slice( viewer, hover_point ) )
		{
			// Delete hovered measurement
			Core::ActionRemove::Dispatch( Core::Interface::GetMouseActionContext(), 
				this->measurements_state_, measurement );

			// Ordering matters here
			this->private_->hover_point_.invalidate();
			this->private_->finish_editing();
			
			return true;
		}
	}

	return false;
}

void MeasurementTool::redraw( size_t viewer_id, const Core::Matrix& proj_mat )
{
	// This function can be called from multiple rendering threads -- one per viewer
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	if ( viewer->is_volume_view() )
	{
		return;
	}
	Core::VolumeSliceHandle vol_slice = viewer->get_active_volume_slice();
	if ( !vol_slice )
	{
		return;
	} 
	
	//-------------- StateEngine locked  -------------------

	// NOTE: The StateEngine and RenderResources should NEVER be locked
	// at the same time because this will lead to deadlock.  So we get all the state we need up
	// front, then unlock the StateEngine, then do the rendering.
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	std::vector< Core::Measurement > measurements = this->measurements_state_->get();
	int active_index = this->active_index_state_->get();
	double opacity = this->opacity_state_->get();
	
	// Getting the length string requires locking the state engine because we need access to the
	// show_world_units state and the LayerManager, which locks the state engine.
	std::vector< std::string > length_strings( measurements.size() );
	for( size_t m_idx = 0; m_idx < measurements.size(); m_idx++ )
	{
		if( measurements[ m_idx ].get_visible() ) 
		{
			length_strings[ m_idx ] = this->get_length_string( measurements[ m_idx ] ); 
		}
		else
		{
			length_strings[ m_idx ] = "";
		}
	}

	size_t active_viewer = static_cast< size_t >( 
		ViewerManager::Instance()->active_viewer_state_->get() );
	
	lock.unlock();
	
	//-------------- StateEngine unlocked  -------------------

	Core::Color in_slice_color = Core::Color( 1.0f, 1.0f, 0.0f );
	Core::Color out_of_slice_color = Core::Color( 0.6f, 0.6f, 0.0f );

	glPushAttrib( GL_LINE_BIT | GL_POINT_BIT | GL_TRANSFORM_BIT );
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd( proj_mat.data() );

	glPointSize( 5.0f );
	glEnable( GL_LINE_SMOOTH );

	// Projected vertices per measurement
	std::vector< std::vector< Core::Point > > vertices( measurements.size() );
	// Are both points in slice for each measurement
	std::vector< bool > both_in_slice( measurements.size() );

	// Draw lines and points
	for( size_t m_idx = 0; m_idx < measurements.size(); m_idx++ )
	{
		if( measurements[ m_idx ].get_visible() )
		{
			Core::Measurement m = measurements[ m_idx ];
			bool vertex_in_slice[ 2 ];

			// TODO: Visually represent point in front differently?
			// Project points, determine if they are "in slice"
			vertices[ m_idx ].resize( 2 );
			for( size_t v_idx = 0; v_idx < 2; v_idx++ )
			{
				Core::Point measurement_point;
				m.get_point( static_cast< int >( v_idx ), measurement_point );

				vertex_in_slice[ v_idx ] = this->private_->in_slice( viewer, measurement_point );
				
				// Project 3D point onto slice
				double x_pos, y_pos;
				vol_slice->project_onto_slice( measurement_point, x_pos, y_pos );
				
				vertices[ m_idx ][ v_idx ][ 0 ] = x_pos;
				vertices[ m_idx ][ v_idx ][ 1 ] = y_pos;
			}

			// Draw line before points so that points are rendered on top (looks better in case 
			// where one point is "in slice", the other is not)
			MeasurementToolPrivate::lock_type lock( this->private_->get_mutex() );
			both_in_slice[ m_idx ] = vertex_in_slice[ 0 ] && vertex_in_slice[ 1 ];
			Core::Color color = both_in_slice[ m_idx ] ? in_slice_color : out_of_slice_color;
			glColor4f( color.r(), color.g(), color.b(), static_cast< float >( opacity ) );

			// If editing this measurement
			if( this->private_->editing_ && 
				m.get_id() == this->private_->hover_point_.measurement_id_ )
			{
				// Draw coarsely dotted line 
				glLineWidth( 2.0f );
				glLineStipple(1, 0x00FF );
				glEnable( GL_LINE_STIPPLE ); 
			}
			else if( static_cast<int>( m_idx ) == active_index )
			{
				// Draw solid line to indicate active measurement
				glLineWidth( 2.0f );
				glDisable( GL_LINE_STIPPLE );
			}
			else
			{
				// Draw finely dotted line
				glLineWidth( 1.5f );
				glLineStipple(1, 0xAAAA );
				glEnable( GL_LINE_STIPPLE ); 
			}

			glBegin( GL_LINES );
			for ( size_t i = 0; i < 2; i++ )
			{
				glVertex2d( vertices[ m_idx ][ i ].x(), vertices[ m_idx ][ i ].y() );
			}   
			glEnd();

			// Draw points
			for( size_t v_idx = 0; v_idx < 2; v_idx++ )
			{
				Core::Color color = vertex_in_slice[ v_idx ] ? in_slice_color : out_of_slice_color;
			
				// Make the moving point red in the current viewer (not necessarily active viewer)
				if( this->private_->viewer_ && 
					viewer_id == this->private_->viewer_->get_viewer_id() && 
					this->private_->editing_ && 
					m.get_id() == this->private_->hover_point_.measurement_id_ &&
					v_idx == this->private_->hover_point_.point_index_ )
				{
					color = Core::Color( 1.0, 0.0, 0.0 );
				}

				// Render GL_POINT
				glColor4f( color.r(), color.g(), color.b(), static_cast< float >( opacity ) );
				glBegin( GL_POINTS );
				glVertex2d( vertices[ m_idx ][ v_idx ].x(), vertices[ m_idx ][ v_idx ].y() );
				glEnd();
			}
		}
	}

	// Render length above line, label below line
	// Different projection matrix is required
	glPopMatrix();

	//-------------- StateEngine locked  -------------------

	// Need to lock state engine, but can't lock it and RenderResources at same time or deadlock
	// could occur.  So do conversion to window coordinates first, then do rendering.
	std::vector< std::vector< Core::Point > > window_vertices( measurements.size() );
	for( size_t m_idx = 0; m_idx < measurements.size(); m_idx++ )
	{
		if( measurements[ m_idx ].get_visible() )
		{
			// NOTE: Y values increase from top to bottom of the window
			
			// Convert vertices to window coords
			Core::Point p0 = vertices[ m_idx ][ 0 ];
			Core::Point p1 = vertices[ m_idx ][ 1 ];
			double p0_x, p0_y, p1_x, p1_y;
			// NOTE: Locks state engine, so has to be called before RenderResources is locked
			viewer->world_to_window( p0.x(), p0.y(), p0_x, p0_y );
			viewer->world_to_window( p1.x(), p1.y(), p1_x, p1_y );
			Core::Point window_p0( p0_x, p0_y, 0.0 );
			Core::Point window_p1( p1_x, p1_y, 0.0 );
			window_vertices[ m_idx ].push_back( window_p0 );
			window_vertices[ m_idx ].push_back( window_p1 );
		}
	}

	//-------------- StateEngine unlocked  -------------------

	//-------------- RenderResources locked  -------------------

	// Dealing with textures, so need to lock RenderResources
	Core::RenderResources::lock_type render_lock( Core::RenderResources::GetMutex() );
	// Create local TextRenderer and Texture2D for each render thread to prevent contention 
	Core::TextRendererHandle text_renderer;
	text_renderer.reset( new Core::TextRenderer );
	Core::Texture2DHandle text_texture;
	text_texture.reset( new Core::Texture2D );
	std::vector< unsigned char > buffer( viewer->get_width() * viewer->get_height(), 0 );

	for( size_t m_idx = 0; m_idx < measurements.size(); m_idx++ )
	{
		if( measurements[ m_idx ].get_visible() )
		{
			Core::Measurement m = measurements[ m_idx ];

			// Find position of label and length

			Core::Point window_p0 = window_vertices[ m_idx ][ 0 ];
			Core::Point window_p1 = window_vertices[ m_idx ][ 1 ];

			// Find mid point
			Core::Point mid_point = window_p0 + ( ( window_p1 - window_p0 ) / 2.0 );

			// Find perpendicular to line, normalize
			Core::Vector measure_normal;
			double dx = window_p1.x() - window_p0.x();
			double dy = window_p0.y() - window_p1.y();
			measure_normal = Core::Vector( dy, dx, 0 );
			measure_normal.normalize();
			int pixel_offset = 10;

			// Put length on one side, label on other
			bool flip_normal = measure_normal.y() < 0;
			Core::Point label_point = mid_point + 
				( measure_normal * pixel_offset * ( flip_normal ? -1 : 1 ) );
			Core::Point length_point = mid_point + 
				( measure_normal * pixel_offset * ( flip_normal ? 1 : -1 ) );

			// Find angle of line
			double angle = 0;
			if( dx == 0 )
			{
				if( dy > 0 )
				{
					angle = 90;
				}
				else
				{
					angle = -90;
				}
			}
			else
			{
				angle = ( atan( dy / dx ) ) * 180 / Core::Pi();
			}
		
			//
			// Label
			//

			int text_height = 20; // Heuristically determined -- computing actual size didn't work
			int text_width = static_cast< int >( m.get_id().size() ) * text_height;

			// Higher for more horizontal lines
			double norm_angle = abs( angle ) / 90.0;
			double x_offset = 0;
			// Label to the left of line
			if( angle < 0 )
			{
				x_offset = norm_angle * -0.5 * text_width; 
			}

			// Shift to the left by some amount so that the label doesn't intersect the line or the 
			// other text
			label_point[ 0 ] = label_point.x() + x_offset;

			double y_offset = 0.5 * text_height;
			label_point[ 1 ] = label_point.y() + y_offset;

			// Render label
			unsigned int font_size = 14; // Matches slice number font size
			text_renderer->render( m.get_id(), &buffer[ 0 ], 
				viewer->get_width(), viewer->get_height(), static_cast< int >( label_point.x() ), 
				viewer->get_height() - static_cast< int >( label_point.y() ), font_size, 0 );

			//
			// Length
			//

			std::string length_string = length_strings[ m_idx ];
			text_width = static_cast< int >( length_string.size() ) * text_height;
			
			x_offset = 0;
			// Label to the left of line
			if( angle > 0 )
			{
				x_offset = norm_angle * -0.5 * text_width; 
			}

			// Shift to the left by some amount so that the label doesn't intersect the line or the 
			// other text
			length_point[ 0 ] = length_point.x() + x_offset;

			// Render length
			text_renderer->render( length_string, &buffer[ 0 ], 
				viewer->get_width(), viewer->get_height(), static_cast< int >( length_point.x() ), 
				viewer->get_height() - static_cast< int >( length_point.y() ), font_size, 0 );
		}
	}

	text_texture->enable();
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	text_texture->set_image( viewer->get_width(), viewer->get_height(),
		GL_ALPHA, &buffer[ 0 ], GL_ALPHA, GL_UNSIGNED_BYTE );

	// Blend the text onto the framebuffer
	glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_ADD );
	glBegin( GL_QUADS );
	glColor4f( in_slice_color.r(), in_slice_color.g(), in_slice_color.b(), 
		static_cast< float >( opacity ) );
	glTexCoord2f( 0.0f, 0.0f );
	glVertex2i( 0, 0 );
	glTexCoord2f( 1.0f, 0.0f );
	glVertex2i( viewer->get_width() - 1, 0 );
	glTexCoord2f( 1.0f, 1.0f );
	glVertex2i( viewer->get_width() - 1, viewer->get_height() - 1 );
	glTexCoord2f( 0.0f, 1.0f );
	glVertex2i( 0, viewer->get_height() - 1 );
	glEnd();
	text_texture->disable();
	
	CORE_CHECK_OPENGL_ERROR();

	glPopAttrib();
	glFinish();

	//-------------- RenderResources unlocked  -------------------
}

bool MeasurementTool::has_2d_visual()
{
	return true;
}

std::string MeasurementTool::get_length_string( const Core::Measurement& measurement ) const
{
	// NOTE: Do not call this function if RendererResources is locked.
	// We need access to the show_world_units_state_, and several functions we call also lock
	// the state engine.  
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	// Thread-safety: We get a handle to the active layer 
	// (get_active_layer is thread-safe), so it can't be deleted out from under 
	// us.  
	LayerHandle active_layer = LayerManager::Instance()->get_active_layer();

	double length = 0;
	bool use_scientific = false;
	
	if( !this->show_world_units_state_->get() && active_layer ) 
	{
		// Index units

		// Convert world units to index units
		// Use grid transform from active layer
		Core::GridTransform grid_transform = active_layer->get_grid_transform();

		// Grid transfrom takes index coords to world coords, so we need inverse
		Core::Transform inverse_transform = grid_transform.get_inverse();

		Core::Point p0, p1;
		measurement.get_point( 0 , p0 );
		measurement.get_point( 1 , p1 );
		Core::Vector measure_vec = p1 - p0;
		measure_vec = inverse_transform.project( measure_vec );
		length = measure_vec.length();

		// Use same formatting policy as status bar for coordinates
		if( 10000 < length ) 
		{
			use_scientific = true;
		}
		else 
		{
			use_scientific = false;
		}
	}
	else
	{
		// World units
		length= measurement.get_length();

		// Use same formatting policy as status bar for coordinates
		if( ( 0.0 < length && length < 0.0001 ) || 1000 < length ) 
		{
			use_scientific = true;
		}
		else 
		{
			use_scientific = false;
		}
	}
	
	// Use same formatting policy as status bar for coordinates
	std::ostringstream oss;
	if( use_scientific ) 
	{
		// Use scientific notation
		oss.precision( 2 );
		oss << std::scientific << length;
		return ( oss.str() );
	}
	else 
	{
		// Format normally
		oss.precision( 3 );
		oss << std::fixed << length;
		return ( oss.str() );
	}
}

} // end namespace Seg3D


