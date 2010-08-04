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

#ifndef APPLICATION_VIEWER_VIEWER_H
#define APPLICATION_VIEWER_VIEWER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// STL includes
#include <map>
#include <vector>

// Boost includes 
#include <boost/smart_ptr.hpp>
#include <boost/function.hpp>
#include <boost/thread.hpp>

// Core includes
#include <Core/Viewer/AbstractViewer.h>
#include <Core/Volume/VolumeSlice.h>
#include <Core/State/State.h>

namespace Seg3D
{

// Forward declarations
class Viewer;
typedef boost::shared_ptr< Viewer > ViewerHandle;
typedef boost::weak_ptr< Viewer > ViewerWeakHandle;

class ViewerPrivate;
typedef boost::shared_ptr< ViewerPrivate > ViewerPrivateHandle;

// Class definition
class Viewer : public Core::AbstractViewer, public boost::enable_shared_from_this< Viewer >
{

	// -- Constructor/Destructor --
public:
	friend class ViewManipulator;
	friend class ViewerPrivate;

	Viewer( size_t viewer_id, bool visible = true, const std::string& mode = Viewer::AXIAL_C );
	virtual ~Viewer();

	// -- Mouse and keyboard events handling --
public:

	typedef boost::function< bool( const Core::MouseHistory&, int, int, int ) > 
		mouse_event_handler_type;
	typedef boost::function< bool( size_t, int, int ) > enter_event_handler_type;
	typedef boost::function< bool( size_t ) > leave_event_handler_type;
	typedef boost::function< bool( int, int, int, int, int ) > 	wheel_event_handler_type;

	virtual void mouse_move_event( const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );
	virtual void mouse_press_event( const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );
	virtual void mouse_release_event( const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );
	virtual void mouse_enter_event( int x, int y );
	virtual void mouse_leave_event();
	virtual bool wheel_event( int delta, int x, int y, int buttons, int modifiers );

	virtual bool key_press_event( int key, int modifiers );

	void set_mouse_move_handler( mouse_event_handler_type func );
	void set_mouse_press_handler( mouse_event_handler_type func );
	void set_mouse_release_handler( mouse_event_handler_type func );
	void set_mouse_enter_handler( enter_event_handler_type func );
	void set_mouse_leave_handler( leave_event_handler_type func );
	void set_wheel_event_handler( wheel_event_handler_type func );
	void reset_mouse_handlers();

	// -- Slice operations --
public:
	// GET_VOLUME_SLICE:
	// Returns the volume slice of the specified layer.
	Core::VolumeSliceHandle get_volume_slice( const std::string& layer_id );

	// GET_ACTIVE_LAYER_SLICE:
	// Returns the volume slice that corresponds to the active layer.
	Core::VolumeSliceHandle get_active_volume_slice() const;

	// MOVE_SLICE_TO:
	// Move the slice to the given world coordinate. Used for picking.
	void move_slice_to( const Core::Point& pt );

private:
	friend class ActionOffsetSlice;

	// OFFSET_SLICE:
	// Offset the slice number by the given value.
	int offset_slice( int delta );

	// -- View information --
public:

	// RESIZE:
	// Set the new size of the viewer.
	virtual void resize( int width, int height );

	// AUTO_VIEW:
	// Auto adjust the view for the active layer
	void auto_view();

	// IS_VOLUME_VIEW:
	// Returns true if the current view mode is volume, otherwise false.
	bool is_volume_view() const;

	// GET_ACTIVE_VIEW_STATE:
	// Returns the view state variable associated with the current view mode.
	Core::StateViewBaseHandle get_active_view_state() const;

	// WINDOW_TO_WORLD:
	// Maps from window coordinates to world coordinates.
	// NOTE: Only call this function when the viewer is in one of the 2D modes.
	void window_to_world( int x, int y, double& world_x, double& world_y ) const;

	// GET_PROJECTION_MATRIX:
	// Get the projection matrix of current view mode.
	// NOTE: Only works in 2D modes.
	void get_projection_matrix( Core::Matrix& proj_mat ) const;

	// UPDATE_STATUS_BAR:
	// Update the status bar to show the data information of the specified layer under
	// the mouse cursor. If no layer is specified, the active layer will be used.
	void update_status_bar( int x, int y, const std::string& layer_id = "" );

	// -- Rendering --
public:
	
	// REDRAW:
	// Emits redraw_signal_.
	void redraw( bool delay_update = false );

	// REDRAW_OVERLAY:
	// Emits redraw_overlay_signal_.
	void redraw_overlay( bool delay_update = false );

	// -- Signals --
public:

	typedef boost::signals2::signal< void( bool ) > redraw_signal_type;

	// REDRAW_SIGNAL:
	// Signals that the scene needs to be redrawn.
	redraw_signal_type redraw_signal_;

	// REDRAW_OVERLAY_SIGNAL_:
	// Signals that the overlay needs to be redrawn.
	redraw_signal_type redraw_overlay_signal_;

	// SLICE_CHANGED_SIGNAL_:
	// Triggered when slice number or viewer visibility is changed.
	// Renderer of other viewers connect to this signal to update the overlay.
	typedef boost::signals2::signal< void ( size_t ) > slice_changed_signal_type;
	slice_changed_signal_type slice_changed_signal_;

	// -- State information --
public:

	Core::StateOptionHandle view_mode_state_;

	Core::StateView2DHandle axial_view_state_;
	Core::StateView2DHandle coronal_view_state_;
	Core::StateView2DHandle sagittal_view_state_;
	Core::StateView3DHandle volume_view_state_;
	Core::StateRangedIntHandle slice_number_state_;

	// 2D viewer state
	Core::StateBoolHandle slice_grid_state_;
	Core::StateBoolHandle slice_visible_state_;
	Core::StateBoolHandle slice_picking_visible_state_;

	// 3D viewer state
	Core::StateBoolHandle volume_slices_visible_state_;
	Core::StateBoolHandle volume_isosurfaces_visible_state_;
	Core::StateBoolHandle volume_volume_rendering_visible_state_;
	Core::StateBoolHandle volume_light_visible_state_;

	Core::StateBoolHandle lock_state_;
	Core::StateBoolHandle overlay_visible_state_;
	Core::StateBoolHandle is_picking_target_state_;

public:
	const static std::string AXIAL_C;
	const static std::string SAGITTAL_C;
	const static std::string CORONAL_C;
	const static std::string VOLUME_C;

private:
	ViewerPrivateHandle private_;

	const static size_t VERSION_NUMBER_C;

};

} // end namespace Seg3D

#endif
