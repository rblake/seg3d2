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

//QtUtils Includes
#include <QtUtils/Bridge/QtBridge.h>

//Qt Gui Includes
#include <Interface/ToolInterface/BinarySmoothDilateErodeFilterInterface.h>
#include "ui_BinarySmoothDilateErodeFilterInterface.h"

//Application Includes
#include <Application/Tools/BinarySmoothDilateErodeFilter.h>

SCI_REGISTER_TOOLINTERFACE( Seg3D, BinarySmoothDilateErodeFilterInterface )

namespace Seg3D
{

class BinarySmoothDilateErodeFilterInterfacePrivate
{
public:
	Ui::BinarySmoothDilateErodeFilterInterface ui_;
};

// constructor
BinarySmoothDilateErodeFilterInterface::BinarySmoothDilateErodeFilterInterface() :
	private_( new BinarySmoothDilateErodeFilterInterfacePrivate )
{
}

// destructor
BinarySmoothDilateErodeFilterInterface::~BinarySmoothDilateErodeFilterInterface()
{
}

// build the interface and connect it to the state manager
bool BinarySmoothDilateErodeFilterInterface::build_widget( QFrame* frame )
{
	// Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );
	this->private_->ui_.horizontalLayout->setAlignment( Qt::AlignHCenter );

	// Step 2 - get a pointer to the tool
	BinarySmoothDilateErodeFilter* tool = dynamic_cast< BinarySmoothDilateErodeFilter* > ( this->tool().get() );
	
	// Step 3 - connect the gui to the tool through the QtBridge
	QtUtils::QtBridge::Connect( this->private_->ui_.target_layer_, 
		tool->target_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.use_active_layer_, 
		tool->use_active_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.erode_, 
		tool->erode_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.dilate_, 
		tool->dilate_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.mask_,
		tool->mask_state_ );	
	QtUtils::QtBridge::Connect( this->private_->ui_.mask_invert_,
		tool->mask_invert_state_ );	

	QtUtils::QtBridge::Connect( this->private_->ui_.replaceCheckBox, 
		tool->replace_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.dilateButton,
		tool->valid_target_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.erodeButton,
		tool->valid_target_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.dilateerodeButton,
		tool->valid_target_state_ );

	QtUtils::QtBridge::Connect( this->private_->ui_.only2d_,
		tool->only2d_state_ );

	QButtonGroup* button_group = new QButtonGroup( this );
	button_group->addButton( this->private_->ui_.active_slice_ );
	button_group->addButton( this->private_->ui_.sagittal_slice_ );
	button_group->addButton( this->private_->ui_.coronal_slice_ );
	button_group->addButton( this->private_->ui_.axial_slice_ );

	QtUtils::QtBridge::Connect( button_group, tool->slice_type_state_ );

	QtUtils::QtBridge::Enable( this->private_->ui_.active_slice_, tool->only2d_state_ );	
	QtUtils::QtBridge::Enable( this->private_->ui_.sagittal_slice_, tool->only2d_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.coronal_slice_, tool->only2d_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.axial_slice_, tool->only2d_state_ );



	QtUtils::QtBridge::Enable( this->private_->ui_.target_layer_, tool->use_active_layer_state_, true );

	QtUtils::QtBridge::Connect( this->private_->ui_.dilateerodeButton, boost::bind( 
		&BinarySmoothDilateErodeFilter::execute_dilateerode, tool, Core::Interface::GetWidgetActionContext() ) );
	QtUtils::QtBridge::Connect( this->private_->ui_.dilateButton, boost::bind( 
		&BinarySmoothDilateErodeFilter::execute_dilate, tool, Core::Interface::GetWidgetActionContext() ) );
	QtUtils::QtBridge::Connect( this->private_->ui_.erodeButton, boost::bind( 
		&BinarySmoothDilateErodeFilter::execute_erode, tool, Core::Interface::GetWidgetActionContext() ) );
	
	return true;

} // end build_widget

} // end namespace Seg3D
