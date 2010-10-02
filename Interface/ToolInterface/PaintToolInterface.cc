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

// boost includes
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <QComboBox>

//Qt Gui Includes
#include <Interface/ToolInterface/PaintToolInterface.h>
#include "ui_PaintToolInterface.h"

//Application Includes
#include <Application/Tools/PaintTool.h>

//QtUtils Includes
#include <QtUtils/Bridge/QtBridge.h>

SCI_REGISTER_TOOLINTERFACE( Seg3D, PaintToolInterface )

namespace Seg3D
{

class PaintToolInterfacePrivate
{
public:
	Ui::PaintToolInterface ui_;

    QtUtils::QtSliderIntCombo *brush_radius_;
	QtUtils::QtSliderDoubleCombo *upper_threshold_;
	QtUtils::QtSliderDoubleCombo *lower_threshold_;
	
};

// constructor
PaintToolInterface::PaintToolInterface() :
	private_( new PaintToolInterfacePrivate )
{
}

// destructor
PaintToolInterface::~PaintToolInterface()
{
}

// build the interface and connect it to the state manager
bool PaintToolInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	this->private_->ui_.setupUi( frame );

	//Add the SliderSpinCombos
	this->private_->brush_radius_ = new QtUtils::QtSliderIntCombo( this, false );
	this->private_->ui_.verticalLayout->addWidget( this->private_->brush_radius_ );

	this->private_->upper_threshold_ = new QtUtils::QtSliderDoubleCombo( this, false );
	this->private_->ui_.upperHLayout_bottom->addWidget( this->private_->upper_threshold_ );
	
	this->private_->lower_threshold_ = new QtUtils::QtSliderDoubleCombo( this, false );
	this->private_->ui_.lowerHLayout_bottom->addWidget( this->private_->lower_threshold_ );
	
	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	PaintTool* tool = dynamic_cast< PaintTool* > ( base_tool_.get() );
	
	//Step 3 - connect the gui to the tool through the QtBridge
	QtUtils::QtBridge::Connect( this->private_->ui_.target_mask_, 
		tool->target_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.use_active_layer_, 
		tool->use_active_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.mask_constraint_, 
		tool->mask_constraint_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.negative_mask_constraint_,
		tool->negative_mask_constraint_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.data_constraint_, 
		tool->data_constraint_layer_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.negative_data_constraint_,
		tool->negative_data_constraint_state_ );
	QtUtils::QtBridge::Connect( this->private_->brush_radius_, 
		tool->brush_radius_state_ );
	QtUtils::QtBridge::Connect( this->private_->upper_threshold_, 
		tool->upper_threshold_state_ );
	QtUtils::QtBridge::Connect( this->private_->lower_threshold_, 
		tool->lower_threshold_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.eraseCheckBox, 
		tool->erase_state_ );
	QtUtils::QtBridge::Connect( this->private_->ui_.show_boundary_,
		tool->show_data_cstr_bound_state_ );
	QtUtils::QtBridge::Enable( this->private_->ui_.target_mask_,
		tool->use_active_layer_state_, true );
	
	boost::function< bool () > condition = boost::lambda::bind( &Core::StateLabeledOption::get, 
		tool->data_constraint_layer_state_.get() ) != Tool::NONE_OPTION_C;
	QtUtils::QtBridge::Enable( this->private_->upper_threshold_, 
		tool->data_constraint_layer_state_, condition );
	QtUtils::QtBridge::Enable( this->private_->lower_threshold_,
		tool->data_constraint_layer_state_, condition );
	
//#if defined ( __APPLE__ )  
//	QFont font;
//	font.setPointSize( 11 );
//	this->setFont( font );
////	this->private_->ui_.target_mask_->setFont( font );
////	this->private_->ui_.mask_constraint_->setFont( font );
////	this->private_->ui_.data_constraint_->setFont( font );
//#endif
	

	return true;
} 

} // end namespace Seg3D

