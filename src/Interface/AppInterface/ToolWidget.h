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

#ifndef INTERFACE_APPINTERFACE_TOOLWIDGET_H
#define INTERFACE_APPINTERFACE_TOOLWIDGET_H

// QT Includes
#include <QtGui/QWidget>
#include <QtGui/QFrame>

// Application includes
#include <Application/Tool/ToolInterface.h>
#include <Application/Tool/ToolFactory.h>

namespace Seg3D
{

class ToolWidget : public QWidget, public ToolInterface
{
Q_OBJECT

// -- constructor/destructor --
public:
	ToolWidget();
	virtual ~ToolWidget();

	// CREATE_WIDGET:
	// The constructor only builds the class. Because this is handled through
	// a factory method we use this auxillary function to build the inner parts
	// of the widget
	bool create_widget( QWidget* parent, ToolHandle& tool );

	// BUILD_WIDGET:
	// Function to create the specific tool widget:
	// This one needs to be overloaded
	virtual bool build_widget( QFrame* frame ) = 0; // << NEEDS TO BE REIMPLEMENTED

	// -- widget internals --
protected:
	QFrame* main_frame_;

};

} //end namespace Seg3D

#endif
