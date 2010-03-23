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
#include <Interface/ToolInterface/IntensityCorrectionFilterInterface.h>
#include "ui_IntensityCorrectionFilterInterface.h"

//Application Includes
#include <Application/Tools/IntensityCorrectionFilter.h>

namespace Seg3D
{

SCI_REGISTER_TOOLINTERFACE(IntensityCorrectionFilterInterface)

class IntensityCorrectionFilterInterfacePrivate
{
public:
	Ui::IntensityCorrectionFilterInterface ui_;
	
    SliderIntCombo *order_;
	SliderDoubleCombo *edge_;
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
	private_->ui_.setupUi( frame );

	//Add the SliderSpinCombos
	this->private_->order_ = new SliderIntCombo();
	private_->ui_.orderHLayout_bottom->addWidget( this->private_->order_ );

	this->private_->edge_ = new SliderDoubleCombo();
	private_->ui_.edgeHLayout_bottom->addWidget( this->private_->edge_ );

	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	IntensityCorrectionFilter* tool =
	    dynamic_cast< IntensityCorrectionFilter* > ( base_tool_.get() );
    
    //Step 3 - set the values for the tool ui from the state engine
	
	    //set default falues for the target option list	
	    std::vector< std::string > temp_option_list = tool->target_layer_state_->option_list();
	    for( size_t i = 0; i < temp_option_list.size(); i++)
	    {   
	        this->private_->ui_.targetComboBox->addItem( QString::fromStdString( temp_option_list[i] ) );
	    } 
        this->private_->ui_.targetComboBox->setCurrentIndex(tool->target_layer_state_->index());
        
        // set the defaults for order
	    int order_min = 0; 
	    int order_max = 0;
	    int order_step = 0;
	    tool->order_state_->get_step( order_step );
	    tool->order_state_->get_range( order_min, order_max );
	    private_->order_->setStep( order_step );
        private_->order_->setRange( order_min, order_max );
        private_->order_->setCurrentValue( tool->order_state_->get() );
        
        // set the defaults for edge
        double edge_min = 0.0; 
	    double edge_max = 0.0;
	    double edge_step = 0.0;
	    tool->edge_state_->get_step( edge_step );
	    tool->edge_state_->get_range( edge_min, edge_max );
	    private_->edge_->setStep( edge_step );
        private_->edge_->setRange( edge_min, edge_max );
        private_->edge_->setCurrentValue( tool->edge_state_->get() ); 
        
        this->private_->ui_.replaceCheckBox->setChecked( tool->replace_state_->get() );


	//Step 4 - connect the gui to the tool through the QtBridge
	QtBridge::Connect( this->private_->ui_.targetComboBox, tool->target_layer_state_ );
	QtBridge::Connect( this->private_->order_, tool->order_state_ );
	QtBridge::Connect( this->private_->edge_, tool->edge_state_ );
	QtBridge::Connect( this->private_->ui_.replaceCheckBox, tool->replace_state_ );

	//Send a message to the log that we have finised with building the tensity Correction Filter Interface
	SCI_LOG_DEBUG("Finished building an Intensity Correction Filter Interface");
	return ( true );

} // end build_widget
} // namespace Seg3D
