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

// QtGui includes
#include "ui_IntensityCorrectionFilterInterface.h"

// Application includes
#include <Application/Tools/IntensityCorrectionFilter.h>

// QtUtils includes
#include <QtUtils/Bridge/QtBridge.h>

// Interface includes
#include <Interface/ToolInterface/IntensityCorrectionFilterInterface.h>


SCI_REGISTER_TOOLINTERFACE( Seg3D, IntensityCorrectionFilterInterface )

namespace Seg3D
{

class IntensityCorrectionFilterInterfacePrivate
{
public:
	Ui::IntensityCorrectionFilterInterface ui_;
};

// constructor
IntensityCorrectionFilterInterface::IntensityCorrectionFilterInterface() :
	private_( new IntensityCorrectionFilterInterfacePrivate )
{
}

// destructor
IntensityCorrectionFilterInterface::~IntensityCorrectionFilterInterface()
{
}

// build the interface and connect it to the state manager
bool IntensityCorrectionFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );
	this->private_->ui_.horizontalLayout_4->setAlignment( Qt::AlignHCenter );
	this->private_->ui_.horizontalLayout_5->setAlignment( Qt::AlignHCenter );

	//Step 2 - get a pointer to the tool
	IntensityCorrectionFilter* tool = dynamic_cast< IntensityCorrectionFilter* > ( 
		this->tool().get() );

	//Step 3 - connect the gui to the tool through the QtBridge
	QtUtils::QtBridge::Connect( this->private_->ui_.target_layer_, 
		tool->target_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.use_active_layer_, 
		tool->use_active_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.preserve_data_format_,
		tool->preserve_data_format_state_ );

	QtUtils::QtBridge::Connect( this->private_->ui_.order_, 
		tool->order_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.edge_, 
		tool->edge_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.replaceCheckBox, 
		tool->replace_state_ );
	
	QtUtils::QtBridge::Enable( this->private_->ui_.runFilterButton, tool->valid_target_state_ );
	QtUtils::QtBridge::Show( this->private_->ui_.message_alert_, tool->valid_target_state_, true );

	QtUtils::QtBridge::Enable( this->private_->ui_.target_layer_, 
		tool->use_active_layer_state_, true );
	QtUtils::QtBridge::Connect( this->private_->ui_.runFilterButton, boost::bind(
		&Tool::execute, tool, Core::Interface::GetWidgetActionContext() ) );
	
	this->private_->ui_.order_->set_description( "Polynomial Order" );
	this->private_->ui_.edge_->set_description( "Edge Sensitivity" );

	return true;
}
	
} // namespace Seg3D
