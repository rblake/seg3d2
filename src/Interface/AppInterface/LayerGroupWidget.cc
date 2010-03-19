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
	
	SliderSpinComboDouble* center_x_adjuster_crop;
	SliderSpinComboDouble* center_y_adjuster_crop;
	SliderSpinComboDouble* center_z_adjuster_crop;
	
	SliderSpinComboDouble* size_width_adjuster_crop;
	SliderSpinComboDouble* size_height_adjuster_crop;
	SliderSpinComboDouble* size_depth_adjuster_crop;
	
    SliderSpinComboDouble* scale_adjuster;
	
	std::string group_id_;
	Utils::GridTransform grid_transform_;
	int current_height;
	int current_width;
	int current_depth;
	//LayerGroupWeakHandle group_handle_;
};
	
LayerGroupWidget::LayerGroupWidget( QWidget* parent, LayerGroupHandle group ) :
	private_( new LayerGroupWidgetPrivate )
{
	this->setParent( parent );
	
	this->private_->ui_.setupUi( this );
	
	this->private_->group_id_ = group->get_group_id();

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
	//this->private_->ui_.open_button_->setChecked( true );
	
	// add the slider spinner combo's for the crop
	this->private_->center_x_adjuster_crop = new SliderSpinComboDouble( this->private_->ui_.widget );
	this->private_->ui_.horizontalLayout_11->addWidget( this->private_->center_x_adjuster_crop );
	this->private_->center_x_adjuster_crop->setObjectName( QString::fromUtf8( "center_x_adjuster_crop" ) );
	
	this->private_->center_y_adjuster_crop = new SliderSpinComboDouble( this->private_->ui_.widget_2 );
	this->private_->ui_.horizontalLayout_12->addWidget( this->private_->center_y_adjuster_crop );
	this->private_->center_y_adjuster_crop->setObjectName( QString::fromUtf8( "center_y_adjuster_crop" ) );
	
	this->private_->center_z_adjuster_crop = new SliderSpinComboDouble( this->private_->ui_.widget_3 );
	this->private_->ui_.horizontalLayout_14->addWidget( this->private_->center_z_adjuster_crop );
	this->private_->center_z_adjuster_crop->setObjectName( QString::fromUtf8( "center_z_adjuster_crop" ) );
	
	this->private_->size_height_adjuster_crop = new SliderSpinComboDouble( this->private_->ui_.widget_4 );
	this->private_->ui_.horizontalLayout_7->addWidget( this->private_->size_height_adjuster_crop );
	this->private_->size_height_adjuster_crop->setObjectName( QString::fromUtf8( "size_height_adjuster_crop" ) );
	
	this->private_->size_width_adjuster_crop = new SliderSpinComboDouble( this->private_->ui_.widget_5 );
	this->private_->ui_.horizontalLayout_9->addWidget( this->private_->size_width_adjuster_crop );
	this->private_->size_width_adjuster_crop->setObjectName( QString::fromUtf8( "size_width_adjuster_crop" ) );
	
	this->private_->size_depth_adjuster_crop = new SliderSpinComboDouble( this->private_->ui_.widget_6 );
	this->private_->ui_.horizontalLayout_10->addWidget( this->private_->size_depth_adjuster_crop );
	this->private_->size_depth_adjuster_crop->setObjectName( QString::fromUtf8( "size_depth_adjuster_crop" ) );
	
	this->private_->scale_adjuster = new SliderSpinComboDouble( this->private_->ui_.widget_7 );
	this->private_->ui_.horizontalLayout_15->addWidget( this->private_->scale_adjuster );
	this->private_->scale_adjuster->setObjectName( QString::fromUtf8( "scale_adjuster" ) );
	
	// set some local values for the current size
	this->private_->current_width = static_cast<int>( group->get_grid_transform().nx() );
	this->private_->current_height = static_cast<int>( group->get_grid_transform().ny() );
	this->private_->current_depth = static_cast<int>( group->get_grid_transform().nz() );
	
	//  connect the gui signals and slots
    connect( this->private_->scale_adjuster, SIGNAL( valueAdjustedContinuously( double ) ), this, SLOT( adjust_new_size_labels( double )) );
	connect( this->private_->ui_.open_button_, SIGNAL( toggled( bool ) ), this, SLOT( show_layers( bool )) );
	connect( this->private_->ui_.group_resample_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_resample( bool )) );
	connect( this->private_->ui_.group_crop_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_crop( bool )) );
	connect( this->private_->ui_.group_transform_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_transform( bool )) );
	connect( this->private_->ui_.group_flip_rotate_button_, SIGNAL( clicked ( bool ) ), this, SLOT( show_flip_rotate( bool )) );
	connect( this->private_->ui_.group_delete_button_, SIGNAL( clicked ( bool ) ), this, SLOT( show_delete( bool )) );
	connect( this->private_->ui_.confirm_delete_checkbox_, SIGNAL( clicked ( bool ) ), this, SLOT( enable_delete_button( bool )) );

	
	// Add all the layers to the group
	layer_list_type temp_layer_list = group->get_layer_list();
	for( layer_list_type::reverse_iterator i = temp_layer_list.rbegin(); i != temp_layer_list.rend(); ++i )
	{
		this->add_layer(( *i ));
	}
	
	//Set the defaulf values for the Group UI and make the connections to the state engine
	    // --- GENERAL ---
	    this->private_->ui_.open_button_->setChecked( group->show_layers_state_.get() );
	    this->private_->ui_.group_visibility_button_->setChecked(group->visibility_state_.get() );
	    
	    QtBridge::Connect( this->private_->ui_.open_button_, group->show_layers_state_ );
	    QtBridge::Connect( this->private_->ui_.group_visibility_button_, group->visibility_state_ );
	
	
	    // --- RESAMPLE ---
	    // = set the default values
	    this->private_->ui_.x_axis_label_current_->setText( QString::fromUtf8("X: ") + QString::number(group->get_grid_transform().nx()) );
	    this->private_->ui_.y_axis_label_current_->setText( QString::fromUtf8("Y: ") + QString::number(group->get_grid_transform().ny()) );
	    this->private_->ui_.z_axis_label_current_->setText( QString::fromUtf8("Z: ") + QString::number(group->get_grid_transform().nz()) );
	    
	    this->private_->ui_.x_axis_label_new_->setText( QString::fromUtf8("X: ") + QString::number(group->get_grid_transform().nx()) );
	    this->private_->ui_.y_axis_label_new_->setText( QString::fromUtf8("Y: ") + QString::number(group->get_grid_transform().ny()) );
	    this->private_->ui_.z_axis_label_new_->setText( QString::fromUtf8("Z: ") + QString::number(group->get_grid_transform().nz()) );
	    
	    QFont font;
        font.setPointSize(8);
        this->private_->ui_.x_axis_label_current_->setFont(font);
        this->private_->ui_.y_axis_label_current_->setFont(font);
        this->private_->ui_.z_axis_label_current_->setFont(font);
        
        font.setBold(true);
        this->private_->ui_.x_axis_label_new_->setFont(font);
        this->private_->ui_.y_axis_label_new_->setFont(font);
        this->private_->ui_.z_axis_label_new_->setFont(font);
        
        this->private_->ui_.resample_replace_checkBox_->setChecked( group->resample_replace_state_->get() );
        
         // = make the connections
        QtBridge::Connect( this->private_->scale_adjuster, group->resample_factor_state_ );
	    QtBridge::Connect( this->private_->ui_.resample_replace_checkBox_, group->resample_replace_state_ );
        
        
        // --- CROP ---
        // = set the default values
        this->private_->size_width_adjuster_crop->setRanges( 0, group->get_grid_transform().nx() );
        this->private_->size_height_adjuster_crop->setRanges( 0, group->get_grid_transform().ny() );
        this->private_->size_depth_adjuster_crop->setRanges( 0, group->get_grid_transform().nz() );
        
        this->private_->size_width_adjuster_crop->setCurrentValue( group->get_grid_transform().nx() );
        this->private_->size_height_adjuster_crop->setCurrentValue( group->get_grid_transform().ny() );
        this->private_->size_depth_adjuster_crop->setCurrentValue( group->get_grid_transform().nz() );

        this->private_->center_x_adjuster_crop->setRanges( 0, group->get_grid_transform().nx() );
        this->private_->center_y_adjuster_crop->setRanges( 0, group->get_grid_transform().ny() );
        this->private_->center_z_adjuster_crop->setRanges( 0, group->get_grid_transform().nz() );
        
        this->private_->center_x_adjuster_crop->setCurrentValue( 0 );
        this->private_->center_y_adjuster_crop->setCurrentValue( 0 );
        this->private_->center_z_adjuster_crop->setCurrentValue( 0 );
        
        this->private_->ui_.crop_replace_checkBox_->setChecked( group->crop_replace_state_->get() );
        
        // = make the connections
        QtBridge::Connect( this->private_->size_width_adjuster_crop, group->crop_size_width_state_ );
        QtBridge::Connect( this->private_->size_height_adjuster_crop, group->crop_size_height_state_ );
        QtBridge::Connect( this->private_->size_depth_adjuster_crop, group->crop_size_depth_state_ );
        
        QtBridge::Connect( this->private_->center_x_adjuster_crop, group->crop_center_x_state_ );
        QtBridge::Connect( this->private_->center_y_adjuster_crop, group->crop_center_y_state_ );
        QtBridge::Connect( this->private_->center_z_adjuster_crop, group->crop_center_z_state_ );
        
        QtBridge::Connect( this->private_->ui_.crop_replace_checkBox_, group->crop_replace_state_ );
        
        
        // --- TRANSFORM ---
        // = set the default values
        this->private_->ui_.spacing_x_spinbox_->setValue( group->get_grid_transform().nx() );
        this->private_->ui_.spacing_y_spinbox_->setValue( group->get_grid_transform().ny() );
        this->private_->ui_.spacing_z_spinbox_->setValue( group->get_grid_transform().nz() );
        
        this->private_->ui_.transform_replace_checkBox_->setChecked( group->resample_replace_state_->get() );
        
        // = make the connections
	    QtBridge::Connect( this->private_->ui_.origin_x_spinbox_, group->transform_origin_x_state_ );
	    QtBridge::Connect( this->private_->ui_.origin_y_spinbox_, group->transform_origin_y_state_ );
	    QtBridge::Connect( this->private_->ui_.origin_z_spinbox_, group->transform_origin_z_state_ );
	    
	    QtBridge::Connect( this->private_->ui_.spacing_x_spinbox_, group->transform_spacing_x_state_ );
	    QtBridge::Connect( this->private_->ui_.spacing_y_spinbox_, group->transform_spacing_y_state_ );
	    QtBridge::Connect( this->private_->ui_.spacing_z_spinbox_, group->transform_spacing_z_state_ );
        
        QtBridge::Connect( this->private_->ui_.transform_replace_checkBox_, group->transform_replace_state_ );
        
        

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

void LayerGroupWidget::adjust_new_size_labels( double scale_factor )
{
    this->private_->ui_.x_axis_label_new_->setText( QString::fromUtf8("X: ") + QString::number(this->private_->current_width * scale_factor ) );
	this->private_->ui_.y_axis_label_new_->setText( QString::fromUtf8("Y: ") + QString::number(this->private_->current_height * scale_factor ) );
	this->private_->ui_.z_axis_label_new_->setText( QString::fromUtf8("Z: ") + QString::number(this->private_->current_depth * scale_factor ) );
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
		this->private_->ui_.group_tools_->show();
	}
	else
	{
		this->private_->ui_.group_frame_->hide();
		this->private_->ui_.group_tools_->hide();
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






