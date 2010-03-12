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

#ifndef INTERFACE_APPINTERFACE_TOOLSDOCKWIDGET_H
#define INTERFACE_APPINTERFACE_TOOLSDOCKWIDGET_H

// QT includes
#include <QtGui>

// STL includes
#include <string>
#include <map>

// Boost includes
#include <boost/signals2/signal.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/shared_ptr.hpp>

// Core includes
#include <Utils/Core/ConnectionHandler.h>

// Application includes
#include <Application/ToolManager/ToolManager.h>

// Interface includes
#include <Interface/AppInterface/ToolWidget.h>
#include <Interface/AppInterface/ToolBoxWidget.h>

namespace Seg3D
{

// Forward declaration
class ToolsDockWidget;

// Class definition
class ToolsDockWidget : public QDockWidget, public Utils::ConnectionHandler
{
Q_OBJECT

// -- constructor/destructor --
public:
	ToolsDockWidget( QWidget *parent = 0 );
	virtual ~ToolsDockWidget();

	// -- functions that control the toolbox --
public:

	// OPEN_TOOL
	// The internals of opening a tool
	void open_tool( ToolHandle& tool );

	// CLOSE_TOOL
	// The internals of closing a tool
	void close_tool( ToolHandle& tool );

	// ACTIVATE_TOOL
	// The internals of activating a tool
	void activate_tool( ToolHandle& tool );

	// -- internals of this class --
private:
	// List of tool widgets
	typedef std::map< std::string, ToolWidget* > tool_widget_list_type;
	tool_widget_list_type tool_widget_list_;

	// Pointer to the ToolBox Widget
	ToolBoxWidget* toolbox_;

	// -- static functions for callbacks into this widget --
public:
	typedef QPointer< ToolsDockWidget > qpointer_type;

	// NOTE: This function do not take references as the parameters are
	// forwarded to a different thread and thus need a copy of the handle
	// to ensure that the tool will not be destructed until this function is
	// handled.

	// HANDLEOPENTOOL:
	// This function should be called to open the tool, this one relays all the
	// information properly to the Qt thread
	static void HandleOpenTool( qpointer_type qpointer, ToolHandle tool );

	// HANDLECLOSETOOL:
	// This function should be called to close the tool, this one relays all the
	// information properly to the Qt thread
	static void HandleCloseTool( qpointer_type qpointer, ToolHandle tool );

	// HANDLEACTIVATETOOL:
	// This function should be called to close the tool, this one relays all the
	// information properly to the Qt thread
	static void HandleActivateTool( qpointer_type qpointer, ToolHandle tool );
};

} // end namespace Seg3D

#endif
