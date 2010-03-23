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
#include <Interface/ToolInterface/HistogramEqualizationFilterInterface.h>
#include "ui_HistogramEqualizationFilterInterface.h"

//Application Includes
#include <Application/Tools/HistogramEqualizationFilter.h>

namespace Seg3D
{

SCI_REGISTER_TOOLINTERFACE(HistogramEqualizationFilterInterface)

class HistogramEqualizationFilterInterfacePrivate
{
public:
	Ui::HistogramEqualizationFilterInterface ui_;
    
    SliderDoubleCombo *upper_threshold_;
	SliderDoubleCombo *lower_threshold_;
	SliderIntCombo *alpha_;
};

// constructor
HistogramEqualizationFilterInterface::HistogramEqualizationFilterInterface() :
	private_( new HistogramEqualizationFilterInterfacePrivate )
{
}

// destructor
HistogramEqualizationFilterInterface::~HistogramEqualizationFilterInterface()
{
}

// build the interface and connect it to the state manager
bool HistogramEqualizationFilterInterface::build_widget( QFrame* frame )
{
	//Step 1 - build the Qt GUI Widget
	private_->ui_.setupUi( frame );

	//Add the SliderSpinCombos
	private_->upper_threshold_ = new SliderDoubleCombo();
	private_->ui_.upperHLayout_bottom->addWidget( private_->upper_threshold_ );

	private_->lower_threshold_ = new SliderDoubleCombo();
	private_->ui_.lowerHLayout_bottom->addWidget( private_->lower_threshold_ );

	private_->alpha_ = new SliderIntCombo();
	private_->ui_.alphaHLayout_bottom->addWidget( private_->alpha_ );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	HistogramEqualizationFilter* tool =
	    dynamic_cast< HistogramEqualizationFilter* > ( base_tool_.get() );
	    
	    
	//Step 3 - set the values for the tool ui from the state engine
	
	    //set default falues for the target option list	
	    std::vector< std::string > temp_option_list = tool->target_layer_state_->option_list();
	    for( size_t i = 0; i < temp_option_list.size(); i++)
	    {   
	        this->private_->ui_.targetComboBox->addItem( QString::fromStdString( temp_option_list[i] ) );
	    } 
        this->private_->ui_.targetComboBox->setCurrentIndex(tool->target_layer_state_->index());
	    
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
        
         // set the defaults for the alpha
	    int alpha_min = 0; 
	    int alpha_max = 0;
	    int alpha_step = 0;
	    tool->alpha_state_->get_step( alpha_step );
	    tool->alpha_state_->get_range( alpha_min, alpha_max );
	    private_->alpha_->setStep( alpha_step );
        private_->alpha_->setRange( alpha_min, alpha_max );
        private_->alpha_->setCurrentValue( tool->alpha_state_->get() );
        
        // set the default for the replace state
        this->private_->ui_.replaceCheckBox->setChecked( tool->replace_state_->get() );
 


	//Step 4 - connect the gui to the tool through the QtBridge
	QtBridge::Connect( private_->ui_.targetComboBox, tool->target_layer_state_ );
	QtBridge::Connect( private_->upper_threshold_, tool->upper_threshold_state_ );
	QtBridge::Connect( private_->lower_threshold_, tool->lower_threshold_state_ );
	QtBridge::Connect( private_->alpha_, tool->alpha_state_ );
	QtBridge::Connect( private_->ui_.replaceCheckBox, tool->replace_state_ );

	//Send a message to the log that we have finised with building the Histogram Equalization Filter Interface
	SCI_LOG_DEBUG("Finished building a Histogram Equalization Filter Interface");
	return ( true );
} // end build_widget

} // end namespace Seg3D
