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
#include <Core/Utils/Log.h>

//QtInterface Includes
#include <QtInterface/Bridge/QtBridge.h>

//Application includes
#include <Application/ProjectManager/ProjectManager.h>

// Interface includes
#include <Interface/AppInterface/ProjectDockWidget.h>
#include "ui_ProjectDockWidget.h"

namespace Seg3D
{

class ProjectDockWidgetPrivate
{
public:

	Ui::ProjectDockWidget ui_;

};

ProjectDockWidget::ProjectDockWidget( QWidget *parent ) :
	QDockWidget( parent ), 
	private_( new ProjectDockWidgetPrivate )
{
	if( this->private_ )
	{
		this->private_->ui_.setupUi( this );
		QtUtils::QtBridge::Connect( this->private_->ui_.project_name_edit_, 
			ProjectManager::Instance()->current_project_->project_name_state_ );
		//this->private_->ui_.project_name_edit_->setText( QString::fromStdString( 
//			ProjectManager::Instance()->current_project_name() ) );
	}
}

ProjectDockWidget::~ProjectDockWidget()
{

}
	


} // end namespace Seg3D
