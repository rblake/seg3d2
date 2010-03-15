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
#include <Interface/ToolInterface/CannyEdgeDetectionFilterInterface.h>
#include "ui_CannyEdgeDetectionFilterInterface.h"

//Application Includes
#include <Application/Tools/CannyEdgeDetectionFilter.h>

namespace Seg3D
{

SCI_REGISTER_TOOLINTERFACE(CannyEdgeDetectionFilterInterface)

class CannyEdgeDetectionFilterInterfacePrivate
{
public:
	Ui::CannyEdgeDetectionFilterInterface ui_;
};

// constructor
CannyEdgeDetectionFilterInterface::CannyEdgeDetectionFilterInterface() :
	private_( new CannyEdgeDetectionFilterInterfacePrivate )
{
}

// destructor
CannyEdgeDetectionFilterInterface::~CannyEdgeDetectionFilterInterface()
{
}

// build the interface and connect it to the state manager
bool CannyEdgeDetectionFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	private_->ui_.setupUi( frame );

	//Add the SliderSpinCombos
	varianceAdjuster = new SliderSpinComboDouble();
	private_->ui_.varianceHLayout_bottom->addWidget( varianceAdjuster );

	errorAdjuster = new SliderSpinComboDouble();
	private_->ui_.errorHLayout_bottom->addWidget( errorAdjuster );

	thresholdAdjuster = new SliderSpinComboDouble();
	private_->ui_.thresholdHLayout_bottom->addWidget( thresholdAdjuster );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	CannyEdgeDetectionFilter* tool = dynamic_cast< CannyEdgeDetectionFilter* > ( base_tool_.get() );

	//Step 3 - connect the gui to the tool through the QtBridge
	QtBridge::Connect( private_->ui_.targetComboBox, tool->target_layer_state_ );
	QtBridge::Connect( varianceAdjuster, tool->variance_state_ );
	QtBridge::Connect( errorAdjuster, tool->max_error_state_ );
	QtBridge::Connect( thresholdAdjuster, tool->threshold_state_ );
	QtBridge::Connect( private_->ui_.replaceCheckBox, tool->replace_state_ );

	//Send a message to the log that we have finised with building the Detection Filter Interface
	SCI_LOG_DEBUG("Finished building a Canny Edge Detection Filter Interface");

	return ( true );
}

} // namespace Seg3D

