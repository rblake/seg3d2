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

// Application includes
#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/MeasurementTool.h>
#include <Application/ViewerManager/ViewerManager.h>

// Core includes
#include <Core/Geometry/Point.h>
#include <Core/State/Actions/ActionAdd.h>
#include <Core/State/Actions/ActionRemove.h>
#include <Core/State/Actions/ActionSetAt.h>
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

	void update_active_index();
	void handle_measurements_changed();
	void handle_units_selection_changed( std::string units );
	void handle_active_layer_changed( LayerHandle active_layer );
	void handle_opacity_changed();
	
	//
	// Viewer 
	//

	// UPDATE_VIEWERS:
	// Redraw in relevant viewers.
	void update_viewers();

	// FIND_HOVER_POINT:
	// Find measurement point (if any) that mouse is currently hovering over.
	bool find_hover_point();

	// MOVE_HOVER_POINT:
	// Move hover point to current mouse position.
	void move_hover_point();

	// UPDATE_HOVER_POINT:
	// Find or move the hover point, depending on whether user is editing or not.
	void update_hover_point();

	bool get_hover_measurement( size_t& index, Core::Measurement& measurement, 
		Core::Point& world_point ); 

	void start_editing();
	void finish_editing();

	void get_mouse_world_point( Core::Point& world_point );

	// IN_SLICE:
	// Return true if the point is in the current slice.  This is slightly different from seed 
	// points -- in this case the point needs to be within epsilon of the slice, not just nearest
	// to it.
	bool in_slice( ViewerHandle viewer, const Core::Point& world_point );

	MeasurementTool* tool_;

	std::string active_group_id_;

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

	int measurement_id_counter_; // Only accessed from interface thread

	Core::MousePosition mouse_pos_; // Only accessed from interface thread
};

void MeasurementToolPrivate::update_active_index()
{
	ASSERT_IS_APPLICATION_THREAD();

	size_t num_measurements = this->tool_->get_measurements().size();
	
	if( num_measurements > 0 )
	{
		// If the active index isn't in the valid range, set to end of list
		int active_index = this->tool_->get_active_index();
		if( active_index == -1 || active_index >= static_cast< int >( num_measurements ) )
		{
			this->tool_->set_active_index( static_cast< int >( num_measurements ) - 1 );
		}
	}
	else
	{
		this->tool_->set_active_index( -1 );
	}
}

void MeasurementToolPrivate::handle_measurements_changed()
{
	// Measurements may have been added or removed, so update the active index
	this->update_active_index();

	// Need to redraw the overlay
	this->update_viewers();
}


void MeasurementToolPrivate::handle_units_selection_changed( std::string units )
{
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
	this->update_viewers();
}

void MeasurementToolPrivate::update_viewers()
{
	/*if ( this->signal_block_count_ > 0 )
	{
		return;
	}*/

	ViewerManager::Instance()->update_2d_viewers_overlay();
}

bool MeasurementToolPrivate::find_hover_point()
{
	lock_type lock( this->get_mutex() );

	// Need to find first point within radius.
	// NOTE: Should we find closest point within radius instead?

	// Compute the size of a pixel in world space
	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	double x0, y0, x1, y1;
	viewer->window_to_world( 0, 0, x0, y0 );
	viewer->window_to_world( 1, 1, x1, y1 );
	double pixel_width = Core::Abs( x1 - x0 );
	double pixel_height = Core::Abs( y1 - y0 );

	// Compute the mouse position in world space
	double world_x, world_y;
	viewer->window_to_world( this->mouse_pos_.x_, this->mouse_pos_.y_, world_x, world_y );

	// Search for the first vertex that's within 2 pixels of current mouse position
	double range_x = pixel_width * 4;
	double range_y = pixel_height * 4;

	Core::VolumeSliceHandle volume_slice = viewer->get_active_volume_slice();

	std::vector< Core::Measurement > measurements = this->tool_->get_measurements();
	BOOST_FOREACH( Core::Measurement m, measurements )
	{
		// Ignore hidden measurements
		if( m.get_visible() )
		{
			// For each measurement point
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

void MeasurementToolPrivate::move_hover_point()
{
	// Should be run from interface thread

	lock_type lock( this->get_mutex() );
	
	if( this->hover_point_.is_valid() )
	{
		// Convert current mouse coords to 3D world point
		Core::Point moved_pt;
		this->get_mouse_world_point( moved_pt );

		// Update hover point
		size_t measurement_index;
		Core::Measurement edited_measurement;
		Core::Point measurement_point;
		// See get_measurement() for thread safety notes
		if( this->get_hover_measurement( measurement_index, edited_measurement, measurement_point ) )
		{
			int point_index = this->hover_point_.point_index_;
			edited_measurement.set_point( point_index, moved_pt );
		}
	
		this->tool_->set_measurement( measurement_index, edited_measurement );
	}
}

void MeasurementToolPrivate::update_hover_point()
{
	lock_type lock( this->get_mutex() );

	if( this->editing_ )
	{
		// Move hovered point.  We don't want to "find" the hover point every time the mouse moves
		// during editing since we might find a different point and start moving that one instead.
		this->move_hover_point();
	}
	else
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();

		// Find the measurement point that the mouse is currently hovering over, if there is one
		if( this->find_hover_point() )
		{
			// Set cross icon color to yellow to indicate that point could be selected (but isn�t)
			// TODO: Change this to something better
			viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );
		}
		else
		{
			// Set cross icon color to grey (normal)
			viewer->set_cursor( Core::CursorShape::CROSS_E );
		}
	}
}

bool MeasurementToolPrivate::get_hover_measurement( size_t& index, Core::Measurement& measurement, 
	Core::Point& world_point )
{
	lock_type lock( this->get_mutex() );

	if( this->hover_point_.is_valid() )
	{
		// Get measurements
		std::vector< Core::Measurement> measurements = this->tool_->get_measurements();

		// Find one with matching id
		for( size_t i = 0; i < measurements.size(); i++ )
		{
			Core::Measurement m = measurements[ i ];
			if( m.get_id() == this->hover_point_.measurement_id_ )
			{
				index = i;
				measurement = m;
				measurement.get_point( this->hover_point_.point_index_, world_point );
				return true;
			}
		}
	}
	return false;
}

void MeasurementToolPrivate::start_editing()
{
	lock_type lock( this->get_mutex() );

	this->editing_ = true;

	// Change cursor to indicate that editing is happening	
	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	// TODO: Change this to something better
	viewer->set_cursor( Core::CursorShape::CLOSED_HAND_E );

	// Redraw so that measurement is drawn with dotted line instead of solid line
	this->update_viewers();
}

void MeasurementToolPrivate::finish_editing()
{
	lock_type lock( this->get_mutex() );

	this->editing_ = false;

	// No need to update measurement (interactively updated) 
	// Change cursor to indicate that editing is not happening
	// NOTE: Cursor has three states: editing, hovering, and normal
	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	
	if( this->hover_point_.is_valid() )
	{
		// TODO: Change this to something better
		viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );
	}
	else
	{
		viewer->set_cursor( Core::CursorShape::CROSS_E );
	}
	
	// Redraw so that measurement is drawn with solid line instead of dotted line
	this->update_viewers();
}

void MeasurementToolPrivate::get_mouse_world_point( Core::Point& world_point )
{
	// Make sure to use current viewer, slice, and mouse coordinates.  These aren't passed
	// as parameters because we don't have this info in the case where the view has been changed
	// by a key command.
	ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
	double world_x, world_y;
	viewer->window_to_world( this->mouse_pos_.x_, this->mouse_pos_.y_, world_x, world_y );
	Core::VolumeSliceHandle active_slice = viewer->get_active_volume_slice();
	active_slice->get_world_coord( world_x, world_y, world_point );
}

bool MeasurementToolPrivate::in_slice( ViewerHandle viewer, const Core::Point& world_point )
{
	// Basically in_slice has to be within epsilon of this slice so that editing won't move the
	// measurement points
	Core::VolumeSliceHandle volume_slice = viewer->get_active_volume_slice();
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
	this->private_->measurement_id_counter_ = 0;

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
	this->add_connection( this->units_selection_state_->value_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_units_selection_changed, this->private_, _2 ) ) );
	this->add_connection( LayerManager::Instance()->active_layer_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_active_layer_changed, this->private_, _1 ) ) );
	this->add_connection( this->opacity_state_->state_changed_signal_.connect(
		boost::bind( &MeasurementToolPrivate::handle_opacity_changed, this->private_ ) ) );
}

MeasurementTool::~MeasurementTool()
{
	this->disconnect_all();
}

std::vector< Core::Measurement > MeasurementTool::get_measurements() const
{
	// NOTE: Need to lock state engine as this function is run from the interface thread
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	return this->measurements_state_->get();

	// TODO:
	// Because indices are used to modify the measurements list, we need to synchronize
	// get and set so that indices don't become invalid in between.
	// Scenario: 
	// - ActionRemove posted by interface thread, but not yet processed on app thread
	// - Interface thread locks state engine, gets measurements list and index of matching point
	// - ActionRemove processed, index no longer valid (or at least incorrect)
	// - ActionSetAt with wrong index
	//Core::ActionGet::Dispatch( Core::Interface::GetWidgetActionContext(), this->measurements_state_ )
}

void MeasurementTool::set_measurement( size_t index, const Core::Measurement& measurement )
{
	// Ensure that state is changed on application thread
	Core::ActionSetAt::Dispatch( Core::Interface::GetMouseActionContext(),
		this->measurements_state_, index, measurement );
}

void MeasurementTool::add_measurement( const Core::Measurement& measurement )
{
	// Add measurement
	// Ensure that state is changed on application thread
	Core::ActionAdd::Dispatch( Core::Interface::GetWidgetActionContext(), 
		this->measurements_state_, measurement );
}

void MeasurementTool::remove_measurement( const Core::Measurement& measurement )
{
	// Set active index to the end of the list 
	size_t num_measurements = this->get_measurements().size();
	if( num_measurements > 1 )
	{
		this->set_active_index( static_cast< int >( num_measurements ) - 1 );
	}

	// Remove measurement
	// Ensure that state is changed on application thread
	Core::ActionRemove::Dispatch( Core::Interface::GetWidgetActionContext(), 
		this->measurements_state_, measurement );
}

int MeasurementTool::get_active_index() const
{
	// NOTE: Need to lock state engine as this function is run from the interface thread
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	return this->active_index_state_->get();
}

void MeasurementTool::set_active_index( int active_index )
{
	// Ensure that state is changed on application thread
	Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(), 
		this->active_index_state_, active_index );
}

bool MeasurementTool::get_show_world_units() const
{
	// NOTE: Need to lock state engine as this function may be run from the interface thread
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	return this->show_world_units_state_->get();
}

double MeasurementTool::get_opacity() const
{
	// NOTE: Need to lock state engine as this function may be run from the rendering thread
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	return this->opacity_state_->get();
}

void MeasurementTool::go_to_active_measurement( int point_index )
{
	// Find 3D point based on active measurement index and point index
	// "Pick" this point
}


bool MeasurementTool::handle_mouse_move( ViewerHandle viewer, 
	const Core::MouseHistory& mouse_history, int button, int buttons, int modifiers )
{
	this->private_->mouse_pos_.x_ = mouse_history.current_.x_;
	this->private_->mouse_pos_.y_ = mouse_history.current_.y_;
	this->private_->update_hover_point();
	
	return false;
}

bool MeasurementTool::handle_mouse_press( ViewerHandle viewer, 
	const Core::MouseHistory& mouse_history, int button, int buttons, int modifiers )
{
	MeasurementToolPrivate::lock_type lock( this->private_->get_mutex() );

	if( button == Core::MouseButton::LEFT_BUTTON_E )
	{
		if( this->private_->editing_ ) 
		{
			// We were editing

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
				// Hovering over point that is in slice 

				// Start editing
				this->private_->start_editing();
			}
			else 
			{
				// Hovering over no point or point not in slice

				// Create new measurement 
				Core::Measurement measurement;

				// Need ID for measurement
				// TODO Find id not already in use (may have been loaded from session file)
				measurement.set_id( "M" + Core::ExportToString( this->private_->measurement_id_counter_ ) );
				this->private_->measurement_id_counter_++;

				// Get 3D point from 2D mouse coords, create new measurement with P1 = P2, add to state vector using action.
				Core::Point pt;
				this->private_->get_mouse_world_point( pt );
				measurement.set_point( 0, pt );
				measurement.set_point( 1, pt );
				measurement.set_visible( true );

				// Add measurement to state vector
				this->add_measurement( measurement );

				// Second point in measurement needs to be the hover point
				this->private_->hover_point_.measurement_id_ = measurement.get_id();
				this->private_->hover_point_.point_index_ = 1; // Editing second point

				this->private_->start_editing();
			}
		}
	}
	//else if( button == Core::MouseButton::MIDDLE_BUTTON_E )
	//{
	//	if( this->private_->hover_point_.is_valid() ) // Hovering
	//	{
	//		// This point does not have to be in current slice
	//		// Snap point to current slice
	//		// TODO Different icon for hovering but not in slice?  Confusing otherwise.
	//	}
	//	else
	//	{
	//		// Default middle-click behavior
	//	}
	//}
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
			this->remove_measurement( measurement );

			// Ordering matters here
			this->private_->hover_point_.invalidate();
			this->private_->finish_editing();
		}
		else
		{
			// Default right-click behavior
			return false;
		}
	}

	return true;
}

bool MeasurementTool::handle_wheel( ViewerHandle viewer, int delta, int x, int y, int buttons, 
	int modifiers )
{	
	// Update hover point based on new slice
	this->private_->update_hover_point();

	return true;
}

bool MeasurementTool::handle_key_press( ViewerHandle viewer, int key, int modifiers )
{
	// If changed view, update hover point based on new view
	this->private_->update_hover_point();

	return true;
}

void MeasurementTool::redraw( size_t viewer_id, const Core::Matrix& proj_mat )
{
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	Core::VolumeSliceHandle vol_slice = viewer->get_active_volume_slice();

	// Apply opacity to color
	double opacity = this->get_opacity();
	Core::Color in_slice_color = Core::Color( 1.0f, 1.0f, 0.0f );
	Core::Color out_of_slice_color = Core::Color( 0.6f, 0.6f, 0.0f );

	glPushAttrib( GL_LINE_BIT | GL_POINT_BIT | GL_TRANSFORM_BIT );
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd( proj_mat.data() );

	glPointSize( 5.0f );
	glLineWidth( 2.0f );
	glEnable( GL_LINE_SMOOTH );
	glLineStipple(1, 0x00FF );

	std::vector< Core::Measurement > measurements = this->get_measurements();
	BOOST_FOREACH( Core::Measurement m, measurements )
	{
		if( m.get_visible() )
		{
			std::vector< Core::Point > vertices( 2 );

			bool both_in_slice = true;

			// TODO Visually represent point in front differently
			// Draw points
			for( size_t i = 0; i < vertices.size(); i++ )
			{
				Core::Point measurement_point;
				m.get_point( static_cast< int >( i ), measurement_point );

				// Set color based on in_slice
				bool in_slice = this->private_->in_slice( viewer, measurement_point );
				Core::Color color = in_slice ? in_slice_color : out_of_slice_color;
			
				if( !in_slice ) 
				{
					both_in_slice = false;
				}
				
				// Project 3D point onto slice
				double x_pos, y_pos;
				vol_slice->project_onto_slice( measurement_point, x_pos, y_pos );
				vertices[ i ][ 0 ] = x_pos;
				vertices[ i ][ 1 ] = y_pos;

				// Render GL_POINT
				glColor4f( color.r(), color.g(), color.b(), static_cast< float >( opacity ) );
				glBegin( GL_POINTS );
				glVertex2d( x_pos, y_pos );
				glEnd();
			}

			// Draw line
			MeasurementToolPrivate::lock_type lock( this->private_->get_mutex() );
			Core::Color color = both_in_slice ? in_slice_color : out_of_slice_color;
			glColor4f( color.r(), color.g(), color.b(), static_cast< float >( opacity ) );

			// If editing this measurement
			if( this->private_->editing_ && 
				m.get_id() == this->private_->hover_point_.measurement_id_ )
			{
				// Draw dotted line (GL_LINE)
				glEnable( GL_LINE_STIPPLE ); 
			}
			else
			{
				// Draw regular line (GL_LINE)
				glDisable( GL_LINE_STIPPLE );
			}

			glBegin( GL_LINES );
			for ( size_t i = 0; i < vertices.size(); i++ )
			{
				glVertex2d( vertices[ i ].x(), vertices[ i ].y() );
			}
			glEnd();

			// Render length above line, label below line -- HOW?
		}
	}
	glPopMatrix();
	glPopAttrib();
}

bool MeasurementTool::has_2d_visual()
{
	return true;
}

} // end namespace Seg3D


