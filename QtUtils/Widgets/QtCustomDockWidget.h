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

#ifndef QTUTILS_WIDGETS_QTCUSTOMDOCKWIDGET_H
#define QTUTILS_WIDGETS_QTCUSTOMDOCKWIDGET_H

// QT includes
#include <QtGui/QDockWidget>
#include <QtGui/QCloseEvent>
#include <QtGui/QKeyEvent>

namespace QtUtils
{
	
class QtCustomDockWidget : public QDockWidget
{
Q_OBJECT
	
Q_SIGNALS:
	void closed();
	
public:
	// - Constructor / Destructor
	QtCustomDockWidget( QWidget *parent = 0 );
	virtual ~QtCustomDockWidget();

	// HIDEEVENT:
	// This function is called by Qt to deliver an event that tells that the
	// widget is being hidden. 
	virtual void closeEvent( QCloseEvent* event );
	
// 	// KEYPRESSEVENT:
// 	// This function is called by Qt when a key is pressed
// 	virtual void keyPressEvent( QKeyEvent* event );
// 	
// 	virtual void mousePressEvent( QMouseEvent * event );
	
};
	
} // end namespace QtUtils

#endif
