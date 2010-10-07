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

//Interface Includes
#include <Interface/ToolInterface/CustomWidgets/TargetComboBox.h>
#include <Interface/ToolInterface/CustomWidgets/MaskComboBox.h>

//Qt Gui Includes
#include <Interface/ToolInterface/ThresholdSegmentationLSFilterInterface.h>
#include "ui_ThresholdSegmentationLSFilterInterface.h"

//Application Includes
#include <Application/Tools/ThresholdSegmentationLSFilter.h>
//#include <Application/Filters/Actions/ActionThresholdSegmentationLS.h>

SCI_REGISTER_TOOLINTERFACE( Seg3D, ThresholdSegmentationLSFilterInterface )

namespace Seg3D
{

class ThresholdSegmentationLSFilterInterfacePrivate
{
public:
	Ui::ThresholdSegmentationLSFilterInterface ui_;
	
// 	QtUtils::QtSliderIntCombo *iterations_;
// 	QtUtils::QtSliderDoubleCombo *upper_threshold_;
// 	QtUtils::QtSliderDoubleCombo *lower_threshold_;
// 	QtUtils::QtSliderDoubleCombo *curvature_;
// 	QtUtils::QtSliderDoubleCombo *edge_;
// 	QtUtils::QtSliderDoubleCombo *propagation_;
// 	TargetComboBox *target_;
// 	MaskComboBox *mask_;
};

// constructor
ThresholdSegmentationLSFilterInterface::ThresholdSegmentationLSFilterInterface() :
	private_( new ThresholdSegmentationLSFilterInterfacePrivate )
{
}

// destructor
ThresholdSegmentationLSFilterInterface::~ThresholdSegmentationLSFilterInterface()
{
}

// build the interface and connect it to the state manager
bool ThresholdSegmentationLSFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );
	this->private_->ui_.horizontalLayout_5->setAlignment( Qt::AlignHCenter );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	ThresholdSegmentationLSFilter* tool =
	    dynamic_cast< ThresholdSegmentationLSFilter* > ( base_tool_.get() );
	    
	//Step 3 - connect the gui to the tool through the QtBridge
	QtUtils::QtBridge::Connect( this->private_->ui_.iterations_, tool->iterations_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.upper_threshold_, tool->upper_threshold_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.lower_threshold_, tool->lower_threshold_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.curvature_, tool->curvature_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.edge_, tool->propagation_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.propagation_, tool->edge_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.replaceCheckBox, tool->replace_state_ );
	
	this->connect( this->private_->ui_.runFilterButton, SIGNAL( clicked() ), this, SLOT( execute_filter() ) );
	
	//Send a message to the log that we have finised with building the Segmentation Level Set Filter Interface
	CORE_LOG_DEBUG("Finished building a Segmentation Level Set Filter Interface");
	return ( true );

} // end build_widget
	
void ThresholdSegmentationLSFilterInterface::enable_run_filter( bool valid )
{
	if( valid )
		this->private_->ui_.runFilterButton->setEnabled( true );
	else
		this->private_->ui_.runFilterButton->setEnabled( false );
}

void ThresholdSegmentationLSFilterInterface::execute_filter()
{
	ToolHandle base_tool_ = tool();
	ThresholdSegmentationLSFilter* tool =
	dynamic_cast< ThresholdSegmentationLSFilter* > ( base_tool_.get() );
	
//	ActionThresholdSegmentationLS::Dispatch( tool->target_layer_state_->export_to_string(), 
//		tool->mask_layer_state_->export_to_string(), tool->iterations_state_->get(),
//		tool->upper_threshold_state_->get(), tool->lower_threshold_state_->get(), 
//		tool->curvature_state_->get(), tool->propagation_state_->get(), tool->edge_state_->get(),
//		tool->replace_state_->get() ); 
}

} // end namespace Seg3D
