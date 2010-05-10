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

// Application includes
#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/Viewer/ViewManipulator.h>

#include <Core/Utils/EnumClass.h>
#include <Core/RendererBase/AbstractRenderer.h>
#include <Core/State/State.h>
#include <Core/Volume/DataVolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

namespace Seg3D
{

// Enums for mouse buttons 
// they have the same values as corresponding Qt ones
SCI_ENUM_CLASS
(
	MouseButton,
	NO_BUTTON_E = 0x00000000,
	LEFT_BUTTON_E = 0x00000001,
	RIGHT_BUTTON_E = 0x00000002,
	MID_BUTTON_E = 0x00000004
)

// Enums for key modifiers
// they have the same values as corresponding Qt ones
SCI_ENUM_CLASS
(
	KeyModifier,
	NO_MODIFIER_E = 0x00000000,
	SHIFT_MODIFIER_E = 0x02000000,
	CONTROL_MODIFIER_E = 0x04000000,
	ALT_MODIFIER_E = 0x08000000
)

class MousePosition
{
public:
	MousePosition( int x_in = 0, int y_in = 0 ) :
		x( x_in ), y( y_in )
	{
	}

	int x;
	int y;
};
typedef boost::shared_ptr< MousePosition > MousePositionHandle;

class MouseHistory
{
public:
	MousePosition left_start;
	MousePosition right_start;
	MousePosition mid_start;
	MousePosition previous;
	MousePosition current;
};

// Forward declarations
class ViewManipulator;
class Viewer;
typedef boost::shared_ptr< Viewer > ViewerHandle;
typedef boost::weak_ptr< Viewer > ViewerWeakHandle;

// Class definition
class Viewer : public Core::StateHandler
{

	// -- constructor/destructor --
public:
	friend class ViewManipulator;
	Viewer( size_t viewer_id );
	virtual ~Viewer();

	// -- mouse events handling --
public:

	typedef boost::function< bool( const MouseHistory&, int, int, int ) > mouse_event_handler_type;

	void mouse_move_event( const MouseHistory& mouse_history, int button, int buttons,
	    int modifiers );
	void mouse_press_event( const MouseHistory& mouse_history, int button, int buttons,
	    int modifiers );
	void mouse_release_event( const MouseHistory& mouse_history, int button, int buttons,
	    int modifiers );
	bool wheel_event( int delta, int x, int y, int buttons, int modifiers );

	void set_mouse_move_handler( mouse_event_handler_type func );
	void set_mouse_press_handler( mouse_event_handler_type func );
	void set_mouse_release_handler( mouse_event_handler_type func );
	void reset_mouse_handlers();

private:
	void update_status_bar( int x, int y );
	void pick_point( int x, int y );
	void adjust_contrast_brightness( int dx, int dy );

	mouse_event_handler_type mouse_move_handler_;
	mouse_event_handler_type mouse_press_handler_;
	mouse_event_handler_type mouse_release_handler_;
	ViewManipulatorHandle view_manipulator_;
	bool adjusting_contrast_brightness_;

public:
	void resize( int width, int height );
	bool is_volume_view() const;
	Core::StateViewBaseHandle get_active_view_state();

protected:
	virtual void state_changed();

	// -- Signals and Slots --
public:

	typedef boost::signals2::signal< void( bool ) > redraw_signal_type;
	redraw_signal_type redraw_signal_;
	redraw_signal_type redraw_overlay_signal_;

	// SLICE_CHANGED_SIGNAL_:
	// Triggered when slice number or viewer visibility is changed.
	// Renderers of other viewers connect to this signal to update the overlay.
	typedef boost::signals2::signal< void ( size_t ) > slice_changed_signal_type;
	slice_changed_signal_type slice_changed_signal_;

	// UPDATE_DISPLAY_SIGNAL_:
	// Triggered when new texture is received.
	typedef boost::signals2::signal< void () > update_display_signal_type;
	update_display_signal_type update_display_signal_;

private:
	void change_view_mode( std::string mode, Core::ActionSource source );
	void set_slice_number( int num, Core::ActionSource source = Core::ActionSource::NONE_E );
	void change_visibility( bool visible, Core::ActionSource source );
	void viewer_lock_state_changed( bool locked );
	void layer_state_changed( int affected_view_modes );

	// -- Data structures for keeping track of slices of layers --
private:
	typedef std::map< std::string, Core::MaskVolumeSliceHandle > mask_slices_map_type;
	typedef std::map< std::string, Core::DataVolumeSliceHandle > data_slices_map_type;

	mask_slices_map_type mask_slices_;
	data_slices_map_type data_slices_;
	Core::VolumeSliceHandle active_layer_slice_;

	void insert_layer( LayerHandle layer );
	void delete_layers( std::vector< LayerHandle > layers );
	void set_active_layer( LayerHandle layer );

public:
	Core::MaskVolumeSliceHandle get_mask_volume_slice( const std::string& layer_id );
	Core::DataVolumeSliceHandle get_data_volume_slice( const std::string& layer_id );

	// -- Mutex and lock --
private:

	typedef boost::recursive_mutex mutex_type;
	typedef boost::unique_lock< mutex_type > lock_type;

	inline mutex_type& get_mutex()
	{
		return this->internal_mutex_;
	}

	mutex_type internal_mutex_;

	// -- Other functions and variables --
public:

	// INSTALL_RENDERER:
	// Install a renderer to the viewer.
	void install_renderer( Core::AbstractRendererHandle renderer );

	// Auto adjust the view for the active layer
	void auto_view();

	// Get the stateid prefix of the viewer
	const std::string& get_stateid() const
	{
		return this->stateid();
	}

	Core::VolumeSliceHandle get_active_layer_slice() const;

	// MOVE_SLICE_TO:
	// Move the slice to the given world coordinate. Used for picking.
	void move_slice_to( const Core::Point& pt );

	// GET_TEXTURE:
	// Returns the texture generated by the renderer
	Core::Texture2DHandle get_texture();

	// GET_OVERLAY_TEXTURE:
	// Returns the overlay texture generated by the renderer
	Core::Texture2DHandle get_overlay_texture();

private:
	
	// Auto adjust the view states so the slices are fully visible
	void adjust_view( Core::VolumeSliceHandle target_slice );

	// Move the active slices to the center of the volume
	void adjust_depth( Core::VolumeSliceHandle target_slice );

	// Auto orient the 3D view for the given slice
	void auto_orient( Core::VolumeSliceHandle target_slice );

	void trigger_redraw( bool delay_update );
	void trigger_redraw_overlay( bool delay_update );
	
	// OFFSET_SLICE:
	// Offset the slice number by the given value.
	void offset_slice( int delta );

	// MOVE_SLICE_BY:
	// Move the active slice by the given offset in world coordinates.
	// Called when the viewer is locked to other viewers of the same mode.
	void move_slice_by( double depth_offset );

	// RESET_ACTIVE_SLICE:
	// Bring the active slice into boundary ( if it's out of boundary ).
	// The active slice can only go out of boundary when the viewer is locked to other viewers.
	void reset_active_slice();

	// SET_TEXTURE:
	// Connected to the redraw_completed_signal_ of the renderer.
	void set_texture( Core::Texture2DHandle texture, bool delay_update );

	// SET_OVERLAY_TEXTURE:
	// Connected to the redraw_overlay_completed_signal_ of the renderer.
	void set_overlay_texture( Core::Texture2DHandle texture, bool delay_update );

private:
	const size_t viewer_id_;
	int width_;
	int height_;

	// SIGNALS_BLOCK_COUNT_:
	// Counts the number of times state change signals are blocked
	// Used when state variables are being changed due to internal program logic.
	size_t signals_block_count_;

	// Counts the number of times the slice number is locked (unchangeable)
	size_t slice_lock_count_;

	typedef std::multimap< std::string, boost::signals2::connection > connection_map_type;
	connection_map_type layer_connection_map_;

	Core::Texture2DHandle texture_;
	Core::Texture2DHandle overlay_texture_;

	// -- State information --
public:

	Core::StateOptionHandle view_mode_state_;

	Core::StateView2DHandle axial_view_state_;
	Core::StateView2DHandle coronal_view_state_;
	Core::StateView2DHandle sagittal_view_state_;
	Core::StateView3DHandle volume_view_state_;

	Core::StateRangedIntHandle slice_number_state_;

	Core::StateBoolHandle slice_grid_state_;
	Core::StateBoolHandle slice_visible_state_;

	Core::StateBoolHandle volume_slices_visible_state_;
	Core::StateBoolHandle volume_isosurfaces_visible_state_;
	Core::StateBoolHandle volume_volume_rendering_visible_state_;

	Core::StateBoolHandle viewer_lock_state_;
	Core::StateBoolHandle viewer_visible_state_;
	Core::StateBoolHandle is_picking_target_state_;

private:
	// Indexed view state variables for quick access
	Core::StateViewBaseHandle view_states_[ 4 ];

public:
	const static std::string AXIAL_C;
	const static std::string SAGITTAL_C;
	const static std::string CORONAL_C;
	const static std::string VOLUME_C;

};

} // end namespace Seg3D

#endif
