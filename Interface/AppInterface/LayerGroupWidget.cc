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

//Boost Includes
#include <boost/lexical_cast.hpp>

//Core Includes - for logging
#include <Utils/Core/Log.h>

//Interface Includes
#include <Interface/QtInterface/QtBridge.h>
#include <Interface/AppInterface/LayerGroupWidget.h>
//#include <Interface/ToolInterface/CustomWidgets/SliderSpinComboInt.h>
//#include <Interface/ToolInterface/CustomWidgets/SliderSpinComboDouble.h>

//UI Includes
#include "ui_LayerGroupWidget.h"

//Application Includes
#include <Application/Layer/DataLayer.h>
#include <Application/Layer/MaskLayer.h>
#include <Application/Layer/LayerGroup.h>


namespace Seg3D
{
	
class LayerGroupWidgetPrivate
{
public:
	Ui::LayerGroupWidget ui_;
	
	SliderSpinComboInt* x_adjuster_crop;
	SliderSpinComboInt* y_adjuster_crop;
	SliderSpinComboInt* z_adjuster_crop;
	
	SliderSpinComboInt* height_adjuster_crop;
	SliderSpinComboInt* width_adjuster_crop;
	SliderSpinComboInt* depth_adjuster_crop;
	
	SliderSpinComboInt* x_adjuster_transform;
	SliderSpinComboInt* y_adjuster_transform;
	SliderSpinComboInt* z_adjuster_transform;
	
	SliderSpinComboDouble* height_adjuster_transform;
	SliderSpinComboDouble* width_adjuster_transform;
	SliderSpinComboDouble* depth_adjuster_transform;
	
	SliderSpinComboInt* scale_adjuster;
	
	std::string group_id_;
	Utils::GridTransform grid_transform_;
	//LayerGroupWeakHandle group_handle_;
};
	
LayerGroupWidget::LayerGroupWidget( QWidget* parent, LayerGroupHandle group ) :
	private_( new LayerGroupWidgetPrivate )
{
	this->setParent( parent );
	
	this->private_->ui_.setupUi( this );

	//this->private_->group_handle_ = group;
	
	// set some values of the GUI

	this->private_->ui_.activate_button_->setText( QString::fromStdString( group->get_grid_transform().get_as_string() ) );

	// hide the tool bars 
	this->private_->ui_.roi_->hide();
	this->private_->ui_.resample_->hide();
	this->private_->ui_.flip_rotate_->hide();
	this->private_->ui_.transform_->hide();
	this->private_->ui_.delete_->hide();
	this->private_->ui_.delete_button_->setEnabled( false );
	this->private_->ui_.open_button_->setChecked( true );
	
	// add the slider spinner combo's for the crop
	this->private_->x_adjuster_crop = new SliderSpinComboInt( this->private_->ui_.widget );
	this->private_->ui_.horizontalLayout_11->addWidget( this->private_->x_adjuster_crop );
	this->private_->x_adjuster_crop->setObjectName( QString::fromUtf8( "x_adjuster_crop" ) );
	
	this->private_->y_adjuster_crop = new SliderSpinComboInt( this->private_->ui_.widget_2 );
	this->private_->ui_.horizontalLayout_12->addWidget( this->private_->y_adjuster_crop );
	this->private_->y_adjuster_crop->setObjectName( QString::fromUtf8( "y_adjuster_crop" ) );
	
	this->private_->z_adjuster_crop = new SliderSpinComboInt( this->private_->ui_.widget_3 );
	this->private_->ui_.horizontalLayout_14->addWidget( this->private_->z_adjuster_crop );
	this->private_->z_adjuster_crop->setObjectName( QString::fromUtf8( "z_adjuster_crop" ) );
	
	this->private_->height_adjuster_crop = new SliderSpinComboInt( this->private_->ui_.widget_4 );
	this->private_->ui_.horizontalLayout_7->addWidget( this->private_->height_adjuster_crop );
	this->private_->height_adjuster_crop->setObjectName( QString::fromUtf8( "height_adjuster_crop" ) );
	
	this->private_->width_adjuster_crop = new SliderSpinComboInt( this->private_->ui_.widget_5 );
	this->private_->ui_.horizontalLayout_9->addWidget( this->private_->width_adjuster_crop );
	this->private_->width_adjuster_crop->setObjectName( QString::fromUtf8( "width_adjuster_crop" ) );
	
	this->private_->depth_adjuster_crop = new SliderSpinComboInt( this->private_->ui_.widget_6 );
	this->private_->ui_.horizontalLayout_10->addWidget( this->private_->depth_adjuster_crop );
	this->private_->depth_adjuster_crop->setObjectName( QString::fromUtf8( "depth_adjuster_crop" ) );
	
	// add the slider spinner combo's for the tranform
	this->private_->x_adjuster_transform = new SliderSpinComboInt( this->private_->ui_.widget_8 );
	this->private_->ui_.horizontalLayout_16->addWidget( this->private_->x_adjuster_transform );
	this->private_->x_adjuster_transform->setObjectName( QString::fromUtf8( "x_adjuster_transform" ) );
	
	this->private_->y_adjuster_transform = new SliderSpinComboInt( this->private_->ui_.widget_9 );
	this->private_->ui_.horizontalLayout_17->addWidget( this->private_->y_adjuster_transform );
	this->private_->y_adjuster_transform->setObjectName( QString::fromUtf8( "y_adjuster_transform" ) );
	
	this->private_->z_adjuster_transform = new SliderSpinComboInt( this->private_->ui_.widget_10 );
	this->private_->ui_.horizontalLayout_18->addWidget( this->private_->z_adjuster_transform );
	this->private_->z_adjuster_transform->setObjectName( QString::fromUtf8( "z_adjuster_transform" ) );
	
	this->private_->height_adjuster_transform = new SliderSpinComboDouble( this->private_->ui_.widget_12 );
	this->private_->ui_.horizontalLayout_20->addWidget( this->private_->height_adjuster_transform );
	this->private_->height_adjuster_transform->setObjectName( QString::fromUtf8( "height_adjuster_transform" ) );
	
	this->private_->width_adjuster_transform = new SliderSpinComboDouble( this->private_->ui_.widget_11 );
	this->private_->ui_.horizontalLayout_19->addWidget( this->private_->width_adjuster_transform );
	this->private_->width_adjuster_transform->setObjectName( QString::fromUtf8( "width_adjuster_transform" ) );
	
	this->private_->depth_adjuster_transform = new SliderSpinComboDouble( this->private_->ui_.widget_13 );
	this->private_->ui_.horizontalLayout_21->addWidget( this->private_->depth_adjuster_transform );
	this->private_->depth_adjuster_transform->setObjectName( QString::fromUtf8( "depth_adjuster_transform" ) );
	
	this->private_->scale_adjuster = new SliderSpinComboInt( this->private_->ui_.widget_7 );
	this->private_->ui_.horizontalLayout_15->addWidget( this->private_->scale_adjuster );
	this->private_->scale_adjuster->setObjectName( QString::fromUtf8( "scale_adjuster" ) );
	
//  connect the signals and slots
	connect( this->private_->ui_.open_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_layers( bool )) );
	connect( this->private_->ui_.group_resample_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_resample( bool )) );
	connect( this->private_->ui_.group_crop_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_crop( bool )) );
	connect( this->private_->ui_.group_transform_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_transform( bool )) );
	connect( this->private_->ui_.group_flip_rotate_button_, SIGNAL( clicked ( bool ) ), this, SLOT( show_flip_rotate( bool )) );
	connect( this->private_->ui_.group_delete_button_, SIGNAL( clicked ( bool ) ), this, SLOT( show_delete( bool )) );
	connect( this->private_->ui_.confirm_delete_checkbox_, SIGNAL( clicked ( bool ) ), this, SLOT( enable_delete_button( bool )) );
//	
	
	// Add all the layers to the group
	layer_list_type temp_layer_list = group->get_layer_list();
	for( layer_list_type::iterator i = temp_layer_list.begin(); i != temp_layer_list.end(); ++i )
	{
		this->add_layer(( *i ));
	}
	
}
	
LayerGroupWidget::~LayerGroupWidget()
{
}
	
void LayerGroupWidget::add_layer( LayerHandle layer )
{
	LayerWidget_handle new_layer_handle( new LayerWidget(this->private_->ui_.group_frame_, layer));
	this->private_->ui_.group_frame_layout_->addWidget( new_layer_handle.data() );
	this->layer_list_.push_back( new_layer_handle );
}

std::string& LayerGroupWidget::get_group_id()
{
	return this->private_->group_id_;
}
	
void LayerGroupWidget::show_selection_checkboxes( bool show )
{
	for( int i = 0; i < layer_list_.size(); ++i)
	{
		layer_list_[i]->show_selection_checkbox( show );
	}
}
	
	
void LayerGroupWidget::enable_delete_button( bool enable )
{
	this->private_->ui_.delete_button_->setEnabled( enable );
}

void LayerGroupWidget::show_layers( bool show )
{
	if (show) 
	{
		this->private_->ui_.group_frame_->show();
	}
	else
	{
		this->private_->ui_.group_frame_->hide();
	}
}

void LayerGroupWidget::show_resample( bool show )
{
	if (show) 
	{
		this->private_->ui_.resample_->show();
		
		this->private_->ui_.roi_->hide();
		this->private_->ui_.group_crop_button_->setChecked( false );
		
		this->private_->ui_.delete_->hide();
		this->private_->ui_.group_delete_button_->setChecked( false );
		
		this->private_->ui_.flip_rotate_->hide();
		this->private_->ui_.group_flip_rotate_button_->setChecked( false );
		
		this->private_->ui_.transform_->hide();
		this->private_->ui_.group_transform_button_->setChecked( false );
		
	}
	else
	{
		this->private_->ui_.resample_->hide();
	}
	show_selection_checkboxes( show );
}

void LayerGroupWidget::show_transform( bool show )
{
	if (show) 
	{
		this->private_->ui_.transform_->show();
		
		this->private_->ui_.roi_->hide();
		this->private_->ui_.group_crop_button_->setChecked( false );
		
		this->private_->ui_.delete_->hide();
		this->private_->ui_.group_delete_button_->setChecked( false );
		
		this->private_->ui_.flip_rotate_->hide();
		this->private_->ui_.group_flip_rotate_button_->setChecked( false );
		
		this->private_->ui_.resample_->hide();
		this->private_->ui_.group_resample_button_->setChecked( false );
	}
	else
	{
		this->private_->ui_.transform_->hide();
	}
	show_selection_checkboxes( show );
}

void LayerGroupWidget::show_crop( bool show )
{
	if (show) 
	{
		this->private_->ui_.roi_->show();
		
		this->private_->ui_.transform_->hide();
		this->private_->ui_.group_transform_button_->setChecked( false );
		
		this->private_->ui_.delete_->hide();
		this->private_->ui_.group_delete_button_->setChecked( false );
		
		this->private_->ui_.flip_rotate_->hide();
		this->private_->ui_.group_flip_rotate_button_->setChecked( false );
		
		this->private_->ui_.resample_->hide();
		this->private_->ui_.group_resample_button_->setChecked( false );
	}
	else
	{
		this->private_->ui_.roi_->hide();
	}
	show_selection_checkboxes( show );
}

void LayerGroupWidget::show_flip_rotate( bool show )
{
	if (show) 
	{
		this->private_->ui_.flip_rotate_->show();
		
		this->private_->ui_.transform_->hide();
		this->private_->ui_.group_transform_button_->setChecked( false );
		
		this->private_->ui_.delete_->hide();
		this->private_->ui_.group_delete_button_->setChecked( false );
		
		this->private_->ui_.roi_->hide();
		this->private_->ui_.group_crop_button_->setChecked( false );
		
		this->private_->ui_.resample_->hide();
		this->private_->ui_.group_resample_button_->setChecked( false );
	}
	else
	{
		this->private_->ui_.flip_rotate_->hide();
	}
	show_selection_checkboxes( show );
}

void LayerGroupWidget::show_delete( bool show )
{
	if (show) 
	{
		this->private_->ui_.delete_->show();
		
		this->private_->ui_.transform_->hide();
		this->private_->ui_.group_transform_button_->setChecked( false );
		
		this->private_->ui_.flip_rotate_->hide();
		this->private_->ui_.group_flip_rotate_button_->setChecked( false );
		
		this->private_->ui_.roi_->hide();
		this->private_->ui_.group_crop_button_->setChecked( false );
		
		this->private_->ui_.resample_->hide();
		this->private_->ui_.group_resample_button_->setChecked( false );
	}
	else
	{
		this->private_->ui_.delete_->hide();
	}
	show_selection_checkboxes( show );
}
	
}  //end namespace Seg3D






