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

// Core includes
#include <Core/Viewer/Mouse.h>
#include <Core/State/StateEngine.h>

// QtUtils includes
#include <QtUtils/Widgets/QtCustomDockWidget.h>

// Application includes
#include <Application/LayerManager/Actions/ActionActivateNextLayer.h>
#include <Application/LayerManager/Actions/ActionActivatePreviousLayer.h>
#include <Application/InterfaceManager/InterfaceManager.h>

namespace QtUtils
{
	
QtCustomDockWidget::QtCustomDockWidget( QWidget *parent ) :
	QDockWidget( parent )
{
}
	
QtCustomDockWidget::~QtCustomDockWidget()
{
}

void QtCustomDockWidget::closeEvent( QCloseEvent* event )
{
	Q_EMIT closed();
	event->accept();
}

void QtCustomDockWidget::keyPressEvent( QKeyEvent* event )
{ 
	int e = 0;
	if( event->key() == Core::Key::KEY_LEFT_E )
	{
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		Seg3D::ActionActivatePreviousLayer::Dispatch( Core::Interface::GetKeyboardActionContext() );
	}
	else if( event->key() == Core::Key::KEY_RIGHT_E )
	{
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		Seg3D::ActionActivateNextLayer::Dispatch( Core::Interface::GetKeyboardActionContext() );
	}
	else
	{
		QWidget::keyPressEvent( event );
	}
}

void QtCustomDockWidget::mousePressEvent( QMouseEvent * event )
{
	this->setFocus();
}


} // end namespace QtUtils
