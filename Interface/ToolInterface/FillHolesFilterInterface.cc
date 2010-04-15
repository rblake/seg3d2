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

//Interface Includes
#include <Interface/QtInterface/QtBridge.h>
#include <Interface/ToolInterface/CustomWidgets/TargetComboBox.h>

//Qt Gui Includes
#include <Interface/ToolInterface/FillHolesFilterInterface.h>
#include "ui_FillHolesFilterInterface.h"

//Application Includes
#include <Application/Tools/FillHolesFilter.h>

namespace Seg3D
{

SCI_REGISTER_TOOLINTERFACE(FillHolesFilterInterface)

class FillHolesFilterInterfacePrivate
{
public:
	Ui::FillHolesFilterInterface ui_;
	TargetComboBox *target_;
};

// constructor
FillHolesFilterInterface::FillHolesFilterInterface() :
	private_( new FillHolesFilterInterfacePrivate )
{
}

// destructor
FillHolesFilterInterface::~FillHolesFilterInterface()
{
}

// build the interface and connect it to the state manager
bool FillHolesFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );
		
		// Add the Combobox
		this->private_->target_ = new TargetComboBox( this );
		this->private_->ui_.activeHLayout->addWidget( this->private_->target_ );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	FillHolesFilter* tool = dynamic_cast< FillHolesFilter* > ( base_tool_.get() );
    
    //Step 3 - set the values for the tool ui from the state engine
	
		// N/A
	

	//Step 4 - connect the gui to the tool through the QtBridge
	QtBridge::Connect( this->private_->target_, tool->target_layer_state_ );

	//Send a message to the log that we have finised with building the Fill Holes Filter Interface"
	SCI_LOG_DEBUG("Finished building a Fill Holes Filter Interface");
	return ( true );
} // end build_widget

} // end namespace Seg3D
