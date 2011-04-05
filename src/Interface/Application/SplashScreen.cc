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

// Boost includes
#include <boost/filesystem.hpp>

// Qt includes
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>

// Core indcludes
#include <Core/Utils/Log.h>

// Interface includes
#include <Interface/Application/SplashScreen.h>
#include <Interface/Application/LayerIOFunctions.h>

// Application includes
#include <Application/ProjectManager/ProjectManager.h>
#include <Application/ProjectManager/Actions/ActionLoadProject.h>
#include <Application/ProjectManager/Actions/ActionNewProject.h>

// Resources includes
#include <Resources/QtResources.h>

#include "ui_SplashScreen.h"

namespace Seg3D
{

class SplashScreenPrivate
{
public:
	Ui::SplashScreen ui_;

};


SplashScreen::SplashScreen( QWidget *parent ) :
	QtUtils::QtCustomDialog( parent ),
	private_( new SplashScreenPrivate )
{
	this->setAttribute( Qt::WA_DeleteOnClose, true );
	
	/*this->setWindowFlags( Qt::Dialog | Qt::CustomizeWindowHint | Qt::WindowTitleHint );*/
	this->setWindowFlags( Qt::Dialog | Qt::FramelessWindowHint );

	this->setModal( true );

	// Set up the private internals of the AppSplash class
	this->private_->ui_.setupUi( this );
	this->setObjectName( QString::fromUtf8( "splashscreen" ) );

	//this->setStyleSheet( QString::fromUtf8( "SplashScreen#splashscreen{ background-color: none; }" ) );

	// Disable these since they arent being used yet.
	this->private_->ui_.load_recent_button_->setEnabled( false );
		
	this->populate_recent_projects();
		
	connect( this->private_->ui_.new_project_button_, SIGNAL( clicked() ), 
		this, SLOT( new_project() ) );
	
	connect( this->private_->ui_.existing_project_button_, SIGNAL( clicked() ), 
		this, SLOT( open_existing() ) );
	
	connect( this->private_->ui_.load_recent_button_, SIGNAL( clicked() ), 
		this, SLOT( open_recent() ) );

	connect( this->private_->ui_.quit_seg3d_button_, SIGNAL( clicked() ), 
		this, SLOT( quit() ) );
	
	connect( this->private_->ui_.quick_open_button_, SIGNAL( clicked() ), 
			this, SLOT( quick_open_file() ) );
	
	connect( this->private_->ui_.recent_project_listwidget_, SIGNAL( itemPressed( QListWidgetItem* ) ),
		this, SLOT( enable_load_recent_button( QListWidgetItem* ) ) );	

	connect( this->private_->ui_.recent_project_listwidget_, 
		SIGNAL( itemDoubleClicked ( QListWidgetItem* ) ),
		this, SLOT( open_recent() ) );	
	
	connect( this, SIGNAL( dialog_closed() ), this->parentWidget(), SLOT( close() ) );

	this->private_->ui_.new_project_button_->setFocus();
}


SplashScreen::~SplashScreen()
{
}

	
void SplashScreen::new_project()
{
	this->new_project_wizard_ = new ProjectWizard( this->parentWidget() );
	connect( this->new_project_wizard_, SIGNAL( canceled() ), this, SLOT( unhide() ) );
	this->new_project_wizard_->show();
	this->hide();
}


void SplashScreen::quit()
{
	reinterpret_cast<QWidget*>( this->parent() )->close();
}


void SplashScreen::unhide()
{
	this->show();
}

	
void SplashScreen::open_existing()
{
	boost::filesystem::path current_projects_path = boost::filesystem::complete( 
		boost::filesystem::path( ProjectManager::Instance()-> get_current_project_folder() ) );


	std::string project_type = std::string( "Open " ) + 
		Core::Application::GetApplicationName()	+ " Project";
		
	std::vector<std::string> project_file_extensions = Project::GetProjectFileExtensions();	
	std::vector<std::string> project_path_extensions = Project::GetProjectPathExtensions();	
	std::string project_file_type =  Core::Application::GetApplicationName() + " Project File (";
	
	for ( size_t j = 0; j < project_file_extensions.size(); j++ )
	{
		project_file_type += std::string( " *" ) + project_file_extensions[ j ];
	}

	for ( size_t j = 0; j < project_path_extensions.size(); j++ )
	{
		project_file_type += std::string( " *" ) + project_path_extensions[ j ];
	}

	project_file_type += " )";

	boost::filesystem::path full_path =  boost::filesystem::path( ( 
		QFileDialog::getOpenFileName ( this, 
		QString::fromStdString( project_type ), 
		QString::fromStdString( current_projects_path.string() ), 
		QString::fromStdString( project_file_type ) ) ).toStdString() ); 


	bool is_path_extension = false;
	for ( size_t j = 0; j < project_path_extensions.size(); j++ )
	{
		if ( boost::filesystem::extension( full_path ) == project_path_extensions[ j ] )
		{
			is_path_extension = true;
			break;
		}
	}


	if ( is_path_extension )
	{
		bool found_s3d_file = false;
		
		if ( boost::filesystem::is_directory( full_path ) )
		{
			boost::filesystem::directory_iterator dir_end;
			for( boost::filesystem::directory_iterator dir_itr( full_path ); 
				dir_itr != dir_end; ++dir_itr )
			{
				std::string filename = dir_itr->filename();
				boost::filesystem::path dir_file = full_path / filename;
				for ( size_t j = 0; j < project_file_extensions.size(); j++ )
				{
					if ( boost::filesystem::extension( dir_file ) ==project_file_extensions[ j ] )
					{
						full_path = dir_file;
						found_s3d_file = true;
						break;
					}
				}
				
				if ( found_s3d_file ) break;
			}
		}
		
		if ( !found_s3d_file )
		{
			QMessageBox::critical( this, 
				"Error reading project file",
				"Error reading project file:\n"
				"The project file is incomplete." );
			return;		
		}
	}
	
	std::string path = full_path.parent_path().string();
	std::string file_name = full_path.filename();

	if( boost::filesystem::exists( full_path ) )
	{
		if ( ! ProjectManager::CheckProjectFile( full_path ) )
		{
			QMessageBox::critical( 0, 
				"Error reading project file",
				"Error reading project file:\n"
				"The project file was saved with newer version of Seg3D" );
			return;
		}		
	
		ActionLoadProject::Dispatch( Core::Interface::GetWidgetActionContext(), path );
		this->close();
	}
}


void SplashScreen::open_recent()
{
	std::vector< RecentProject > recent_projects;
	ProjectManager::Instance()->get_recent_projects_from_database( recent_projects );

	int index = this->private_->ui_.recent_project_listwidget_->currentRow();

	if( recent_projects[ index ].name_ == 
		( this->private_->ui_.recent_project_listwidget_->currentItem()->text() ).toStdString() )
	{
		std::vector<std::string> file_project_extensions = Project::GetProjectFileExtensions(); 
	
		bool found_project_file = false;
		boost::filesystem::path path;
		
		for ( size_t j = 0; j < file_project_extensions.size(); j++ )
		{
			path = boost::filesystem::path( recent_projects[ index ].path_ ) / 
			( recent_projects[ index ].name_ + file_project_extensions[ j ] );
			if ( boost::filesystem::exists( path ) ) 
			{
				found_project_file = true;
				break;
			}
		}
	
		if ( ! found_project_file ) return;
		
		if ( ! ProjectManager::CheckProjectFile( path ) )
		{
			QMessageBox::critical( 0, 
				"Error reading project file",
				"Error reading project file:\n"
				"The project file was saved with newer version of Seg3D" );
			return;
		}

		ActionLoadProject::Dispatch( Core::Interface::GetWidgetActionContext(),
			path.string() );
	
	this->close();
	
	}
}


void SplashScreen::quick_open_file()
{
	ActionNewProject::Dispatch( Core::Interface::GetWidgetActionContext(), "",
		"Untitled Project" );
	this->hide();
	LayerIOFunctions::ImportFiles( dynamic_cast< QMainWindow* >( this->parentWidget() ), "" );
	
}


void SplashScreen::populate_recent_projects()
{
	std::vector< RecentProject > recent_projects;
	ProjectManager::Instance()->get_recent_projects_from_database( recent_projects );

	for( size_t i = 0; i < recent_projects.size(); ++i )
	{
		if( recent_projects[ i ].name_ != "" )
		{
			QListWidgetItem *new_item;
			new_item = new QListWidgetItem( QString::fromStdString( recent_projects[ i ].name_ ) );
			new_item->setToolTip( QString::fromUtf8( "This project was created on: " ) +
				QString::fromStdString( recent_projects[ i ].date_ )
				+ QString::fromUtf8( " and is located at: " ) 
				+ QString::fromStdString( recent_projects[ i ].path_ ) );
			this->private_->ui_.recent_project_listwidget_->addItem( new_item );

		}
	}
		
}

	
void SplashScreen::enable_load_recent_button( QListWidgetItem* not_used )
{
	this->private_->ui_.load_recent_button_->setEnabled( true );
}
	
} // end namespace Seg3D
