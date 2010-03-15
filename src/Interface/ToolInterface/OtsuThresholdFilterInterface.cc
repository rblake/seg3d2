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

//Qt Gui Includes
#include <Interface/ToolInterface/OtsuThresholdFilterInterface.h>
#include "ui_OtsuThresholdFilterInterface.h"

//Application Includes
#include <Application/Tools/OtsuThresholdFilter.h>

namespace Seg3D
{

SCI_REGISTER_TOOLINTERFACE(OtsuThresholdFilterInterface)

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
}

// build the interface and connect it to the state manager
bool OtsuThresholdFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	private_->ui_.setupUi( frame );

	// add sliderspincombos
	orderAdjuster = new SliderSpinComboInt();
	private_->ui_.orderHLayout_bottom->addWidget( orderAdjuster );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	OtsuThresholdFilter* tool = dynamic_cast< OtsuThresholdFilter* > ( base_tool_.get() );

	//Step 3 - connect the gui to the tool through the QtBridge
	QtBridge::Connect( private_->ui_.targetComboBox, tool->target_layer_state_ );
	QtBridge::Connect( orderAdjuster, tool->order_state_ );

	//Send a message to the log that we have finised with building the Otsu Threshold Filter Interface
	SCI_LOG_DEBUG("Finished building an Otsu Threshold Filter Interface");
	return ( true );

} // end build_widget

} // end namespace Seg3D
