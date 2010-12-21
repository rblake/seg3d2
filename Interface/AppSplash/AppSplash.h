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

#ifndef INTERFACE_APPSPLASH_APPSPLASH_H
#define INTERFACE_APPSPLASH_APPSPLASH_H

// QT includes
#include <QtCore/QPointer>
#include <QtGui/QListWidgetItem>

// Boost includes
#include <boost/shared_ptr.hpp>

// QtUtils includes
#include <QtUtils/Widgets/QtCustomDialog.h>

// Interface includes
#include <Interface/AppProjectWizard/AppProjectWizard.h>

namespace Seg3D
{

class AppSplashPrivate;

class AppSplash : public QtUtils::QtCustomDialog
{
Q_OBJECT
	
Q_SIGNALS:
	void dialog_closed();

public:
	AppSplash( QWidget* parent = 0 );
	virtual ~AppSplash();
	
private:
	
	// POPULATE_RECENT_PROJECTS:
	// this function is called in the constructor to populate the list of recent projects from the
	// project manager
	void populate_recent_projects();

private Q_SLOTS:
	
	// NEW_PROJECT:
	// this function is called when a user clicks the button to create a new project.  It starts the
	// new project wizard and closes the splash screen
	void new_project();
	
	// OPEN_EXISTING:
	// this function is called when a user clicks the button to open an existing project.  It starts
	// a file browse dialog.
	void open_existing();
	
	// OPEN_RECENT:
	// this function is called when a user double-clicks on an item in the list. 
	void open_recent();

	// quit:
	// this function close the application
	void quit();

	// QUICK_OPEN_FILE:
	// this function is called when someone clicks the quick open file button	
	void quick_open_file();
	
	// ENABLE_LOAD_RECENT_BUTTON:
	// this function enables the load recent button if you select a recent project from the list
	void enable_load_recent_button( QListWidgetItem* item );
	void unhide();

private:
	// Internals of the splash screen
	boost::shared_ptr< AppSplashPrivate > private_;
	
	std::vector< std::string > recent_project_list_;
	QPointer< AppProjectWizard > new_project_wizard_;
	

};

} // end namespace

#endif // APPSPLASH_H
