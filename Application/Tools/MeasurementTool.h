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

#ifndef APPLICATION_TOOLS_MEASUREMENTTOOL_H
#define APPLICATION_TOOLS_MEASUREMENTTOOL_H

// Application includes
#include <Application/Tool/Tool.h>


// Core includes
#include <Core/Geometry/Measurement.h>
#include <Core/State/StateVector.h>

namespace Seg3D
{

class MeasurementTool;
typedef boost::shared_ptr< MeasurementTool > MeasurementToolHandle;

class MeasurementToolPrivate;
typedef boost::shared_ptr< MeasurementToolPrivate > MeasurementToolPrivateHandle;

class MeasurementTool : public Tool
{
friend class MeasurementTableModel;

SEG3D_TOOL
(
	SEG3D_TOOL_NAME( "MeasurementTool", "Tool for creating measurements in slices" )
	SEG3D_TOOL_MENULABEL( "Measure --TEST--" )
	SEG3D_TOOL_MENU( "Tools" )
	SEG3D_TOOL_SHORTCUT_KEY( "Ctrl+ALT+0" )
	SEG3D_TOOL_URL( "http://www.sci.utah.edu/SCIRunDocs/index.php/CIBC:Seg3D2:MeasurementTool:1" )
)

	// -- constructor/destructor --
public:
	MeasurementTool( const std::string& toolid );
	virtual ~MeasurementTool();

	// HANDLE_MOUSE_MOVE:
	// Called when the mouse moves in a viewer.
	virtual bool handle_mouse_move( ViewerHandle viewer, 
		const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );

	// HANDLE_MOUSE_PRESS:
	// Called when a mouse button has been pressed.
	virtual bool handle_mouse_press( ViewerHandle viewer, 
		const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );

	// REDRAW:
	// Draw seed points in the specified viewer.
	// The function should only be called by the renderer, which has a valid GL context.
	virtual void redraw( size_t viewer_id, const Core::Matrix& proj_mat );
	virtual bool has_2d_visual();

	// GET_LENGTH_STRING:
	// Get length as formatted string, respecting index vs. world units.  Function is here so that
	// both interface and rendering use consistent formatting of length.
	std::string get_length_string( const Core::Measurement& measurement ) const;
	
	// -- signals --
public:
	// UNITS_CHANGED_SIGNAL:
	typedef boost::signals2::signal< void () > units_changed_signal_type;
	units_changed_signal_type units_changed_signal_;

	// NUM_MEASUREMENTS_CHANGED_SIGNAL:
	typedef boost::signals2::signal< void () > num_measurements_changed_signal_type;
	units_changed_signal_type num_measurements_changed_signal_;

	// -- state --
public:
	// List of measurements
	Core::StateMeasurementVectorHandle measurements_state_;

	// Index of active measurement in table
	Core::StateIntHandle active_index_state_;

	// Selection between display of index and world units.  Needed for radio button group.
	Core::StateLabeledOptionHandle units_selection_state_; 

	// Boolean indicating whether world units (true) or index units (false) should be displayed.
	Core::StateBoolHandle show_world_units_state_; 
	
	// The opacity of all the measurements for this tool
	Core::StateRangedDoubleHandle opacity_state_;

public:
	static const std::string INDEX_UNITS_C;
	static const std::string WORLD_UNITS_C;

private:
	MeasurementToolPrivateHandle private_;
};

} // end namespace

#endif
