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
#include <Interface/ToolInterface/PaintToolInterface.h>
#include "ui_PaintToolInterface.h"
//#include <Interface/ToolInterface/CustomWidgets/SliderIntCombo.h>
//#include <Interface/ToolInterface/CustomWidgets/SliderDoubleCombo.h>
//#include <Interface/ToolInterface/CustomWidgets/SliderSpinComboInt.h>
//#include <Interface/ToolInterface/CustomWidgets/SliderSpinComboDouble.h>

//Application Includes
#include <Application/Tools/PaintTool.h>

namespace Seg3D
{

SCI_REGISTER_TOOLINTERFACE(PaintToolInterface)

class PaintToolInterfacePrivate
{
public:
	Ui::PaintToolInterface ui_;

    SliderIntCombo *brush_radius_;
	SliderDoubleCombo *upper_threshold_;
	SliderDoubleCombo *lower_threshold_;
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
	private_->ui_.setupUi( frame );

	    //Add the SliderSpinCombos
	    private_->brush_radius_ = new SliderIntCombo( this, true );
	    private_->ui_.verticalLayout->addWidget( private_->brush_radius_ );

        private_->upper_threshold_ = new SliderDoubleCombo();
	    private_->ui_.upperHLayout_bottom->addWidget( private_->upper_threshold_ );
	    
	    private_->lower_threshold_ = new SliderDoubleCombo();
	    private_->ui_.lowerHLayout_bottom->addWidget( private_->lower_threshold_ );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	PaintTool* tool = dynamic_cast< PaintTool* > ( base_tool_.get() );
	
	//Step 3 - set the values for the tool ui from the state engine
	
	    //set default falues for the target option list	
	    std::vector< std::string > temp_option_list = tool->target_layer_state_->option_list();
	    for( size_t i = 0; i < temp_option_list.size(); i++)
	    {   
	        this->private_->ui_.targetComboBox->addItem( QString::fromStdString( temp_option_list[i] ) );
	    } 
        this->private_->ui_.targetComboBox->setCurrentIndex(tool->target_layer_state_->index());
	    
	    //set default falues for the target option list	
	    temp_option_list = tool->mask_layer_state_->option_list();
	    for( size_t i = 0; i < temp_option_list.size(); i++)
	    {   
	        this->private_->ui_.maskComboBox->addItem( QString::fromStdString( temp_option_list[i] ) );
	    } 
        this->private_->ui_.maskComboBox->setCurrentIndex(tool->target_layer_state_->index());
	    
	    
	    // set the defaults for the paint brush size
	    int brush_min = 0; 
	    int brush_max = 0;
	    int brush_radius_step = 0;
	    tool->brush_radius_state_->get_step( brush_radius_step );
	    tool->brush_radius_state_->get_range( brush_min, brush_max );
	    private_->brush_radius_->setStep( brush_radius_step );
        private_->brush_radius_->setRange( brush_min, brush_max );
        private_->brush_radius_->setCurrentValue( tool->brush_radius_state_->get() );
        
        // set the defaults for the upper threshold
        double upper_threshold_min = 0.0; 
	    double upper_threshold_max = 0.0;
	    double upper_threshold_step = 0.0;
	    tool->upper_threshold_state_->get_step( upper_threshold_step );
	    tool->upper_threshold_state_->get_range( upper_threshold_min, upper_threshold_max );
	    private_->upper_threshold_->setStep( upper_threshold_step );
        private_->upper_threshold_->setRange( upper_threshold_min, upper_threshold_max );
        private_->upper_threshold_->setCurrentValue( tool->upper_threshold_state_->get() );
        
        // set the defaults for the lower threshold
        double lower_threshold_min = 0.0; 
	    double lower_threshold_max = 0.0;
	    double lower_threshold_step = 0.0;
	    tool->lower_threshold_state_->get_step( lower_threshold_step );
	    tool->lower_threshold_state_->get_range( lower_threshold_min, lower_threshold_max );
	    private_->lower_threshold_->setStep( lower_threshold_step );
        private_->lower_threshold_->setRange( lower_threshold_min, lower_threshold_max );
        private_->lower_threshold_->setCurrentValue( tool->lower_threshold_state_->get() );
        
        // set the default setchecked state
        private_->ui_.eraseCheckBox->setChecked( tool->erase_state_->get() );
	
	//Step 4 - connect the gui to the tool through the QtBridge
	QtBridge::Connect( private_->ui_.targetComboBox, tool->target_layer_state_ );
	QtBridge::Connect( private_->ui_.maskComboBox, tool->mask_layer_state_ );
	QtBridge::Connect( private_->brush_radius_, tool->brush_radius_state_ );
	QtBridge::Connect( private_->upper_threshold_, tool->upper_threshold_state_ );
	QtBridge::Connect( private_->lower_threshold_, tool->lower_threshold_state_ );
	QtBridge::Connect( private_->ui_.eraseCheckBox, tool->erase_state_ );
	
	
   	//Send a message to the log that we have finised with building the Paint Brush Interface
	SCI_LOG_MESSAGE("Finished building a Paint Brush Interface");

	return ( true );
} // end build_widget

} // end namespace Seg3D

