/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2014 Scientific Computing and Imaging Institute,
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
#include <Core/Interface/Interface.h>

// QtGui includes
#include "ui_RBFFilterInterface.h"

//Application Includes
#include <rbf/Application/Tools/RBFFilter.h>

// QtUtils includes
#include <QtUtils/Bridge/QtBridge.h>

//Interface Includes
#include <rbf/Interface/ToolInterface/RBFFilterInterface.h>


SCI_REGISTER_TOOLINTERFACE( Seg3D, RBFFilterInterface )

namespace Seg3D
{

class RBFFilterInterfacePrivate
{
public:
	Ui::RBFFilterInterface ui_;	
};

// constructor
RBFFilterInterface::RBFFilterInterface() :
	private_( new RBFFilterInterfacePrivate )
{
}

// destructor
RBFFilterInterface::~RBFFilterInterface()
{
}

// build the interface and connect it to the state manager
bool RBFFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );

	//Step 2 - get a pointer to the tool
	RBFFilter* tool = dynamic_cast< RBFFilter* > ( this->tool().get() );
	    
	//Step 3 - connect the gui to the tool through the QtBridge
	QtUtils::QtBridge::Connect( this->private_->ui_.targetLayer_, tool->target_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.useActiveLayer_, tool->use_active_layer_state_ );
	
//	QtUtils::QtBridge::Connect( this->private_->ui_.isovalue_, tool->isovalue_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.sampleX_, tool->sample_x_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.sampleY_, tool->sample_y_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.sampleZ_, tool->sample_z_state_ );
//	QtUtils::QtBridge::Connect( this->private_->ui_.cap_isosurface_, tool->cap_isosurface_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.kernel_, tool->kernel_state_ );

	QtUtils::QtBridge::Connect( this->private_->ui_.runFilterButton_, boost::bind(
		&Tool::execute, tool, Core::Interface::GetWidgetActionContext() ) );
	QtUtils::QtBridge::Connect( this->private_->ui_.clearSeedsButton_, boost::bind(
   &SeedPointsTool::clear, tool, Core::Interface::GetWidgetActionContext() ) );

	QtUtils::QtBridge::Show( this->private_->ui_.messageAlert_, tool->valid_target_state_, true );

	QtUtils::QtBridge::Enable( this->private_->ui_.runFilterButton_, tool->valid_target_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.targetLayer_, tool->use_active_layer_state_, true );

//	this->private_->ui_.isovalue_->set_description( "Isovalue" );
	this->private_->ui_.sampleX_->set_description( "Sample in X" );
	this->private_->ui_.sampleY_->set_description( "Sample in Y" );
	this->private_->ui_.sampleZ_->set_description( "Sample in Z" );


  // TODO: tmp
//	QtUtils::QtBridge::Enable( this->private_->ui_.isovalue_, tool->isovalue_state_, false );
//	QtUtils::QtBridge::Enable( this->private_->ui_.cap_isosurface_, tool->cap_isosurface_state_, false );

//	boost::function< bool () > condition = boost::lambda::bind( &Core::StateLabeledOption::get,
//    tool->targetLayer_state_.get() ) != Tool::NONE_OPTION_C;
//	QtUtils::QtBridge::Enable( this->private_->ui_.isovalue_, tool->targetLayer_state_, condition );

	return true;
}
	
} // end namespace Seg3D
