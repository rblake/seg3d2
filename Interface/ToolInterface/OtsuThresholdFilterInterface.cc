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

// Qt includes
#include <QPointer>

// QtGui includes
#include "ui_OtsuThresholdFilterInterface.h"

// Application includes
#include <Application/Tools/OtsuThresholdFilter.h>
#include <Application/LayerManager/LayerManager.h>

// QtUtils includes
#include <QtUtils/Bridge/QtBridge.h>

// Interface includes
#include <Interface/ToolInterface/OtsuThresholdFilterInterface.h>

// Core includes
#include <Core/DataBlock/Histogram.h>

SCI_REGISTER_TOOLINTERFACE( Seg3D, OtsuThresholdFilterInterface )

namespace Seg3D
{

class OtsuThresholdFilterInterfacePrivate
{
public:
	Ui::OtsuThresholdFilterInterface ui_;

};

// constructor
OtsuThresholdFilterInterface::OtsuThresholdFilterInterface() :
	private_( new OtsuThresholdFilterInterfacePrivate )
{
}

// destructor
OtsuThresholdFilterInterface::~OtsuThresholdFilterInterface()
{
	this->disconnect_all();
}

// build the interface and connect it to the state manager
bool OtsuThresholdFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );

	this->private_->ui_.horizontalLayout_3->setAlignment( Qt::AlignHCenter );

	//Step 2 - get a pointer to the tool
	OtsuThresholdFilter* tool = dynamic_cast< OtsuThresholdFilter* > ( this->tool().get() );
	
	// Step 3 - Qt connections
	{
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );	
		this->private_->ui_.target_layer_->setDisabled( tool->use_active_layer_state_->get() );

		this->connect( this->private_->ui_.use_active_layer_, SIGNAL( toggled( bool ) ),
			this->private_->ui_.target_layer_, SLOT( setDisabled( bool ) ) );

		this->connect( this->private_->ui_.runFilterButton, SIGNAL( clicked() ), 
			this, SLOT( run_filter() ) );

		this->connect( this->private_->ui_.target_layer_, SIGNAL( currentIndexChanged( QString ) ), 
			this, SLOT( refresh_histogram( QString ) ) );
	}

	//Step 4 - connect the gui to the tool through the QtBridge
	QtUtils::QtBridge::Connect( this->private_->ui_.target_layer_, 
		tool->target_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.use_active_layer_, 
		tool->use_active_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.amount_, 
		tool->amount_state_ );
	
	QtUtils::QtBridge::Enable( this->private_->ui_.runFilterButton,
		tool->valid_target_state_ );
	
	QtUtils::QtBridge::Show( this->private_->ui_.message_alert_, tool->valid_target_state_, true );
		
	QtUtils::QtBridge::Enable( this->private_->ui_.histogram_, tool->valid_target_state_ );

	return true;
	
} // end build_widget
	
void OtsuThresholdFilterInterface::run_filter()
{
	tool()->execute( Core::Interface::GetWidgetActionContext() );
}

void OtsuThresholdFilterInterface::refresh_histogram( QString layer_name )
{
	if( layer_name == "" || 
		layer_name == Tool::NONE_OPTION_C.c_str() )
	{
		return;
	}

	DataLayerHandle data_layer = boost::dynamic_pointer_cast< DataLayer >(
		LayerManager::Instance()->get_layer_by_name( layer_name.toStdString() ) );
	if ( !data_layer )
	{
		return;
	}

	this->private_->ui_.histogram_->set_histogram( data_layer->get_data_volume()->
		get_data_block()->get_histogram() );	
}


} // end namespace Seg3D
