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

#ifndef INTERFACE_TOOLINTERFACE_MEASUREMENTTOOLINTERFACE_H
#define INTERFACE_TOOLINTERFACE_MEASUREMENTTOOLINTERFACE_H

// Qt includes
#include <QtCore/QPointer>

// Boost includes
#include <boost/shared_ptr.hpp>

// Base class of the tool widget include
#include <Interface/Application/ToolWidget.h>

namespace Seg3D
{

class MeasurementToolInterfacePrivate;
typedef boost::shared_ptr< MeasurementToolInterfacePrivate > MeasurementToolInterfacePrivateHandle;

// Forward declaration
class MeasurementToolInterface;

class MeasurementToolInterface : public ToolWidget
{
Q_OBJECT


public:
	// Constructor/destructor
	MeasurementToolInterface();
	virtual ~MeasurementToolInterface();
	virtual bool build_widget( QFrame* frame );

private Q_SLOTS:
	// HANDLE_LINEEDIT_NAME_CHANGED:
	// Set the active measurement name state based on the line edit.
	void handle_name_lineedit_changed();

	// HANDLE_TEXTBOX_COMMENT_CHANGED:
	void handle_comment_textbox_changed();

	// HANDLE_LINEEDIT_LENGTH_CHANGED:
	void handle_length_lineedit_changed();

	// HANDLE_CHECKBOX_SHOW_CHANGED:
	void handle_show_checkbox_changed();

	// HANDLE_COLOR_BUTTON_CHANGED:
	void handle_color_button_clicked();

private:
    MeasurementToolInterfacePrivateHandle private_;

public:
	typedef QPointer< MeasurementToolInterface > qpointer_type;

	// UPDATEMEASUREMENTTABLE:
	// Update entire table including dimensions.  Scroll to active index.  
	// Slower than UpdateMeasurementCells, so use only when needed. De-selects selected rows.  
	// TODO: This may be slow due to resizeColumns[Rows]ToContents -- try hard-coding sizes.
	static void UpdateGeneralTab( qpointer_type measurement_interface );

	// UPDATEMEASUREMENTMODEL:
	// Update only table cells, not table dimensions.  Does not scroll to active index.
	static void UpdateTableCells( qpointer_type measurement_interface );

	// UPDATEACTIVEINDEX:
	// Update interface table and text box in response to changed active index.
	// Locks: StateEngine
	static void UpdateTableActiveIndex( qpointer_type measurement_interface );

	// UPDATEACTIVEMEASUREMENTTAB:
	// Set widget values based on current measurement state.
	static void UpdateActiveTab( qpointer_type measurement_interface );
};

} // end namespace Seg3D

#endif
