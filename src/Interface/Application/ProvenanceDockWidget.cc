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

// STL includes
#include <time.h>

// Core includes
#include <Core/Utils/Log.h>
#include <Core/Interface/Interface.h>
#include <Core/State/Actions/ActionSet.h>

// Qt includes
#include <QtGui/QStandardItemModel>
#include <QtGui/QMessageBox>

//QtUtils Includes
#include <QtUtils/Bridge/QtBridge.h>
#include <QtUtils/Utils/QtPointer.h>

//Application includes
#include <Application/ProjectManager/ProjectManager.h>


// Interface includes
#include <Interface/Application/ProvenanceDockWidget.h>
#include "ui_ProvenanceDockWidget.h"

namespace Seg3D
{

class ProvenanceDockWidgetPrivate
{
public:
	Ui::ProvenanceDockWidget ui_;
	bool resetting_;
	std::vector< std::string > sessions_;

};

ProvenanceDockWidget::ProvenanceDockWidget( QWidget *parent ) :
	QtUtils::QtCustomDockWidget( parent ), 
	private_( new ProvenanceDockWidgetPrivate )
{
	if( this->private_ )
	{
		this->private_->ui_.setupUi( this );
		
		this->add_connection( ProjectManager::Instance()->current_project_changed_signal_.
			connect( boost::bind( &ProvenanceDockWidget::HandleProjectChanged, qpointer_type( this ) ) ) );
		
		this->connect_project();
	}

}

ProvenanceDockWidget::~ProvenanceDockWidget()
{
	this->disconnect_all();
}
	

	
void ProvenanceDockWidget::HandleProvenanceResult( qpointer_type qpointer, 
	std::vector< std::pair< ProvenanceID, std::string > > provenance_trail )
{

	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind( 
		&ProvenanceDockWidget::populate_provenance_list, qpointer.data(), provenance_trail ) ) );

}
	
void ProvenanceDockWidget::HandleProjectChanged( qpointer_type qpointer )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind( 
		&ProvenanceDockWidget::connect_project, qpointer.data() ) ) );
	
}
	

void ProvenanceDockWidget::populate_provenance_list( 
	std::vector< std::pair< ProvenanceID, std::string > > provenance_trail )
{
	if( this->isHidden() ) this->show();
	this->raise();

	this->private_->ui_.provenance_list_->verticalHeader()->setUpdatesEnabled( false );

	// Clean out the old table
	for( int j = this->private_->ui_.provenance_list_->rowCount() -1; j >= 0; j-- )
	{
		this->private_->ui_.provenance_list_->removeRow( j );
	}
	
	for( size_t i = 0; i < provenance_trail.size(); ++i )
	{
		QTableWidgetItem *new_prov_number = new QTableWidgetItem( QString::number( provenance_trail[ i ].first ) );
		QTableWidgetItem *new_prov_action = new QTableWidgetItem( QString::fromStdString( provenance_trail[ i ].second ) );

		
		this->private_->ui_.provenance_list_->insertRow( static_cast< int >( i ) );
		this->private_->ui_.provenance_list_->setItem( static_cast< int >( i ), 0, new_prov_number );
		this->private_->ui_.provenance_list_->setItem( static_cast< int >( i ), 1, new_prov_action );
		this->private_->ui_.provenance_list_->verticalHeader()->resizeSection( static_cast< int >( i ), 24 );
	}

	this->private_->ui_.provenance_list_->verticalHeader()->setUpdatesEnabled( true );
	this->private_->ui_.provenance_list_->repaint();
}
	
void ProvenanceDockWidget::connect_project()
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	// Grab this one, once again, but now within the lock
	ProjectHandle current_project = ProjectManager::Instance()->get_current_project();
	if( !current_project ) return;
	
	add_connection( current_project->provenance_record_signal_.
		connect( boost::bind( &ProvenanceDockWidget::HandleProvenanceResult, 
		qpointer_type( this ), _1 ) ) );
}	


} // end namespace Seg3D
