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
#include <Interface/AppInterface/StyleSheet.h>

//UI Includes
#include "ui_LayerGroupWidget.h"

//Application Includes
#include <Application/Layer/DataLayer.h>
#include <Application/Layer/MaskLayer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/LayerManager/Actions/ActionDeleteLayers.h>
#include <Application/LayerManager/Actions/ActionNewMaskLayer.h>



namespace Seg3D
{
	
class LayerGroupWidgetPrivate
{
public:
	Ui::LayerGroupWidget ui_;
	
	SliderDoubleCombo* center_x_adjuster_crop_;
	SliderDoubleCombo* center_y_adjuster_crop_;
	SliderDoubleCombo* center_z_adjuster_crop_;
	
	SliderDoubleCombo* size_width_adjuster_crop_;
	SliderDoubleCombo* size_height_adjuster_crop_;
	SliderDoubleCombo* size_depth_adjuster_crop_;
	
    SliderDoubleCombo* scale_adjuster_;
};
	
LayerGroupWidget::LayerGroupWidget( QWidget* parent, LayerHandle layer ) :
	QWidget( parent ),
	private_( new LayerGroupWidgetPrivate )
{	
    LayerGroupHandle group = layer->get_layer_group();

/*	this->setParent( parent );*/
	
	this->private_->ui_.setupUi( this );
	
	this->group_id_ = group->get_group_id();
	
	// Set up the Drag and Drop
	this->setAcceptDrops( true );

	// set some values of the GUI
	std::string group_name = Utils::ToString( group->get_grid_transform().get_nx() ) + " x " +
		Utils::ToString( group->get_grid_transform().get_ny() ) + " x " +
		Utils::ToString( group->get_grid_transform().get_nz() );
	this->private_->ui_.activate_button_->setText( QString::fromStdString( group_name ) );

	// hide the tool bars 
	this->private_->ui_.roi_->hide();
	this->private_->ui_.resample_->hide();
	this->private_->ui_.flip_rotate_->hide();
	this->private_->ui_.transform_->hide();
	this->private_->ui_.delete_->hide();
	this->private_->ui_.delete_button_->setEnabled( false );
	
	// add the slider spinner combo's for the crop
	this->private_->center_x_adjuster_crop_ = new SliderDoubleCombo( this->private_->ui_.widget );
	this->private_->ui_.horizontalLayout_11->addWidget( this->private_->center_x_adjuster_crop_ );
	this->private_->center_x_adjuster_crop_->setObjectName( QString::fromUtf8( "center_x_adjuster_crop_" ) );
	
	this->private_->center_y_adjuster_crop_ = new SliderDoubleCombo( this->private_->ui_.widget_2 );
	this->private_->ui_.horizontalLayout_12->addWidget( this->private_->center_y_adjuster_crop_ );
	this->private_->center_y_adjuster_crop_->setObjectName( QString::fromUtf8( "center_y_adjuster_crop_" ) );
	
	this->private_->center_z_adjuster_crop_ = new SliderDoubleCombo( this->private_->ui_.widget_3 );
	this->private_->ui_.horizontalLayout_14->addWidget( this->private_->center_z_adjuster_crop_ );
	this->private_->center_z_adjuster_crop_->setObjectName( QString::fromUtf8( "center_z_adjuster_crop_" ) );
	
	this->private_->size_height_adjuster_crop_ = new SliderDoubleCombo( this->private_->ui_.widget_4 );
	this->private_->ui_.horizontalLayout_7->addWidget( this->private_->size_height_adjuster_crop_ );
	this->private_->size_height_adjuster_crop_->setObjectName( QString::fromUtf8( "size_height_adjuster_crop_" ) );
	
	this->private_->size_width_adjuster_crop_ = new SliderDoubleCombo( this->private_->ui_.widget_5 );
	this->private_->ui_.horizontalLayout_9->addWidget( this->private_->size_width_adjuster_crop_ );
	this->private_->size_width_adjuster_crop_->setObjectName( QString::fromUtf8( "size_width_adjuster_crop_" ) );
	
	this->private_->size_depth_adjuster_crop_ = new SliderDoubleCombo( this->private_->ui_.widget_6 );
	this->private_->ui_.horizontalLayout_10->addWidget( this->private_->size_depth_adjuster_crop_ );
	this->private_->size_depth_adjuster_crop_->setObjectName( QString::fromUtf8( "size_depth_adjuster_crop_" ) );
	
	this->private_->scale_adjuster_ = new SliderDoubleCombo( this->private_->ui_.widget_7 );
	this->private_->ui_.horizontalLayout_15->addWidget( this->private_->scale_adjuster_ );
	this->private_->scale_adjuster_->setObjectName( QString::fromUtf8( "scale_adjuster_" ) );
	
	// set some local values for the current size
	this->current_width_ = static_cast<int>( group->get_grid_transform().get_nx() );
	this->current_height_ = static_cast<int>( group->get_grid_transform().get_ny() );
	this->current_depth_ = static_cast<int>( group->get_grid_transform().get_nz() );

	//  connect the gui signals and slots
    connect( this->private_->scale_adjuster_, SIGNAL( valueAdjusted( double ) ), this, SLOT( adjust_new_size_labels( double )) );
	connect( this->private_->ui_.open_button_, SIGNAL( toggled( bool ) ), this, SLOT( show_layers( bool )) );
	connect( this->private_->ui_.group_resample_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_resample( bool )) );
	connect( this->private_->ui_.group_crop_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_crop( bool )) );
	connect( this->private_->ui_.group_transform_button_, SIGNAL( clicked( bool ) ), this, SLOT( show_transform( bool )) );
	connect( this->private_->ui_.group_flip_rotate_button_, SIGNAL( clicked ( bool ) ), this, SLOT( show_flip_rotate( bool )) );
	connect( this->private_->ui_.group_delete_button_, SIGNAL( clicked ( bool ) ), this, SLOT( show_delete( bool )) );
	connect( this->private_->ui_.confirm_delete_checkbox_, SIGNAL( clicked ( bool ) ), this, SLOT( enable_delete_button( bool )) );
	connect( this->private_->ui_.delete_button_, SIGNAL( clicked () ), this, SLOT( uncheck_delete_confirm() ) );
	
	// Add the current layer to the new group
	this->insert_layer( layer, -1 );

	
	//Set the defaulf values for the Group UI and make the connections to the state engine
	    // --- GENERAL ---
	    this->private_->ui_.open_button_->setChecked( group->show_layers_state_.get() );
	    this->private_->ui_.group_visibility_button_->setChecked( group->visibility_state_.get() );
	    
	    QtBridge::Connect( this->private_->ui_.open_button_, group->show_layers_state_ );
	    QtBridge::Connect( this->private_->ui_.group_visibility_button_, group->visibility_state_ );
	    QtBridge::Connect( this->private_->ui_.delete_button_, boost::bind( &ActionDeleteLayers::Dispatch, group ) );
	    QtBridge::Connect( this->private_->ui_.group_new_button_, boost::bind( &ActionNewMaskLayer::Dispatch, group ) );
	
	
	    // --- RESAMPLE ---
	    // = set the default values
	    this->private_->ui_.x_axis_label_current_->setText( QString::fromUtf8("X: ") + 
			QString::number( group->get_grid_transform().get_nx()) );
	    this->private_->ui_.y_axis_label_current_->setText( QString::fromUtf8("Y: ") + 
			QString::number( group->get_grid_transform().get_ny()) );
	    this->private_->ui_.z_axis_label_current_->setText( QString::fromUtf8("Z: ") + 
			QString::number( group->get_grid_transform().get_nz()) );
	    
	    this->private_->ui_.x_axis_label_new_->setText( QString::fromUtf8("X: ") + 
			QString::number( group->get_grid_transform().get_nx()) );
	    this->private_->ui_.y_axis_label_new_->setText( QString::fromUtf8("Y: ") + 
			QString::number( group->get_grid_transform().get_ny()) );
	    this->private_->ui_.z_axis_label_new_->setText( QString::fromUtf8("Z: ") + 
			QString::number( group->get_grid_transform().get_nz()) );
	    
	    QFont font;
        font.setPointSize( 10 );
        this->private_->ui_.x_axis_label_current_->setFont( font );
        this->private_->ui_.y_axis_label_current_->setFont( font );
        this->private_->ui_.z_axis_label_current_->setFont( font );
        
        font.setBold( true );
        this->private_->ui_.x_axis_label_new_->setFont( font );
        this->private_->ui_.y_axis_label_new_->setFont( font );
        this->private_->ui_.z_axis_label_new_->setFont( font );
        
        this->private_->ui_.resample_replace_checkBox_->setChecked( 
			group->resample_replace_state_->get() );
        
        // set the defaults for the upper threshold
        double resample_min = 0.0; 
	    double resample_max = 0.0;
	    double resample_step = 0.0;
	    group->resample_factor_state_->get_step( resample_step );
	    group->resample_factor_state_->get_range( resample_min, resample_max );
	    this->private_->scale_adjuster_->setStep( resample_step );
        this->private_->scale_adjuster_->setRange( resample_min, resample_max );
        this->private_->scale_adjuster_->setCurrentValue( group->resample_factor_state_->get() );
        
         // = make the connections
        QtBridge::Connect( this->private_->scale_adjuster_, group->resample_factor_state_ );
	    QtBridge::Connect( this->private_->ui_.resample_replace_checkBox_, 
			group->resample_replace_state_ );
        

        // --- CROP ---
        // = set the default values
        this->private_->size_width_adjuster_crop_->setRange( 0, 
			group->get_grid_transform().get_nx() );
        this->private_->size_height_adjuster_crop_->setRange( 0, 
			group->get_grid_transform().get_ny() );
        this->private_->size_depth_adjuster_crop_->setRange( 0, 
			group->get_grid_transform().get_nz() );
        
        this->private_->size_width_adjuster_crop_->setCurrentValue( 
			group->get_grid_transform().get_nx() );
        this->private_->size_height_adjuster_crop_->setCurrentValue( 
			group->get_grid_transform().get_ny() );
        this->private_->size_depth_adjuster_crop_->setCurrentValue( 
			group->get_grid_transform().get_nz() );

        this->private_->center_x_adjuster_crop_->setRange( 0, group->get_grid_transform().get_nx() );
        this->private_->center_y_adjuster_crop_->setRange( 0, group->get_grid_transform().get_ny() );
        this->private_->center_z_adjuster_crop_->setRange( 0, group->get_grid_transform().get_nz() );
        
        this->private_->center_x_adjuster_crop_->setCurrentValue( 0 );
        this->private_->center_y_adjuster_crop_->setCurrentValue( 0 );
        this->private_->center_z_adjuster_crop_->setCurrentValue( 0 );
        
        this->private_->ui_.crop_replace_checkBox_->setChecked( group->crop_replace_state_->get() );
        
        // = make the connections
        QtBridge::Connect( this->private_->size_width_adjuster_crop_, group->crop_size_width_state_ );
        QtBridge::Connect( this->private_->size_height_adjuster_crop_, group->crop_size_height_state_ );
        QtBridge::Connect( this->private_->size_depth_adjuster_crop_, group->crop_size_depth_state_ );
        
        QtBridge::Connect( this->private_->center_x_adjuster_crop_, group->crop_center_x_state_ );
        QtBridge::Connect( this->private_->center_y_adjuster_crop_, group->crop_center_y_state_ );
        QtBridge::Connect( this->private_->center_z_adjuster_crop_, group->crop_center_z_state_ );
        
        QtBridge::Connect( this->private_->ui_.crop_replace_checkBox_, group->crop_replace_state_ );
        
        
        // --- TRANSFORM ---
        // = set the default values
        this->private_->ui_.spacing_x_spinbox_->setValue( group->get_grid_transform().get_nx() );
        this->private_->ui_.spacing_y_spinbox_->setValue( group->get_grid_transform().get_ny() );
        this->private_->ui_.spacing_z_spinbox_->setValue( group->get_grid_transform().get_nz() );
        
        this->private_->ui_.transform_replace_checkBox_->setChecked( group->resample_replace_state_->get() );
        
        // = make the connections
	    QtBridge::Connect( this->private_->ui_.origin_x_spinbox_, group->transform_origin_x_state_ );
	    QtBridge::Connect( this->private_->ui_.origin_y_spinbox_, group->transform_origin_y_state_ );
	    QtBridge::Connect( this->private_->ui_.origin_z_spinbox_, group->transform_origin_z_state_ );
	    
	    QtBridge::Connect( this->private_->ui_.spacing_x_spinbox_, group->transform_spacing_x_state_ );
	    QtBridge::Connect( this->private_->ui_.spacing_y_spinbox_, group->transform_spacing_y_state_ );
	    QtBridge::Connect( this->private_->ui_.spacing_z_spinbox_, group->transform_spacing_z_state_ );
        
        QtBridge::Connect( this->private_->ui_.transform_replace_checkBox_, group->transform_replace_state_ );

	this->private_->ui_.group_frame_layout_->setAlignment( Qt::AlignTop );

}
	
LayerGroupWidget::~LayerGroupWidget()
{
}
	
void LayerGroupWidget::insert_layer( LayerHandle layer, int index )
{
	LayerWidgetQHandle new_layer_handle( new LayerWidget(this->private_->ui_.group_frame_, layer ) );
	
	if( index == -1 )
	{
		// If the layer is a data layer, put it on the bottom.
		if( new_layer_handle->get_volume_type() == Utils::VolumeType::DATA_E )
			this->private_->ui_.group_frame_layout_->insertWidget( -1, new_layer_handle.data() );
		else
			this->private_->ui_.group_frame_layout_->insertWidget( 0, new_layer_handle.data() );
	}
	else
	{
		this->private_->ui_.group_frame_layout_->insertWidget( index, new_layer_handle.data() );
	}
	this->layer_list_.push_back( new_layer_handle );
	this->repaint();
}

	
bool LayerGroupWidget::delete_layer( LayerHandle layer )
{	
	for( QVector< LayerWidgetQHandle >::iterator i = this->layer_list_.begin(); i != 
		this->layer_list_.end(); ++i)
	{
		if( ( *i )->get_layer_id() == layer->get_layer_id() )
		{
			( *i )->deleteLater();
			this->layer_list_.erase( i );
			this->repaint();
			return true;
		}
	}
	return false;  
}
	

LayerWidgetQWeakHandle LayerGroupWidget::set_active_layer( LayerHandle layer )
{
	std::string layer_id = layer->get_layer_id();
    for( int i = 0; i < this->layer_list_.size(); ++i)
	{
	    if( layer_id == this->layer_list_[i]->get_layer_id() )
	    {
	        this->layer_list_[i]->set_active( true );
	        this->set_active( true );
			return this->layer_list_[i];
	    }
	}
	return LayerWidgetQWeakHandle();
	
}

void  LayerGroupWidget::set_active( bool active )
{
	if( active )
    {
        this->private_->ui_.base_->setStyleSheet( 
			StyleSheet::GROUP_WIDGET_BASE_ACTIVE_C );
		            	
		this->private_->ui_.group_background_->setStyleSheet( 
			StyleSheet::GROUP_WIDGET_BACKGROUND_ACTIVE_C );
		                
		this->private_->ui_.activate_button_->setStyleSheet( 
			StyleSheet::GROUP_WIDGET_ACTIVATE_BUTTON_ACTIVE_C);	
    }
    else
    {
        this->private_->ui_.base_->setStyleSheet( 
			StyleSheet::GROUP_WIDGET_BASE_INACTIVE_C );
        
        this->private_->ui_.group_background_->setStyleSheet( 
			StyleSheet::GROUP_WIDGET_BACKGROUND_INACTIVE_C );
		              
	    this->private_->ui_.activate_button_->setStyleSheet( 
			StyleSheet::GROUP_WIDGET_ACTIVATE_BUTTON_INACTIVE_C );               
    }
}


const std::string& LayerGroupWidget::get_group_id()
{
	return this->group_id_;
}

	
void LayerGroupWidget::show_selection_checkboxes( bool show )
{
	for( int i = 0; i < this->layer_list_.size(); ++i)
	{
		this->layer_list_[i]->show_selection_checkbox( show );
		this->layer_list_[i]->set_group_menu_status( show );
		this->layer_list_[i]->setAcceptDrops( !show ); 
	}
}

void LayerGroupWidget::adjust_new_size_labels( double scale_factor )
{
    this->private_->ui_.x_axis_label_new_->setText( QString::fromUtf8("X: ") + 
		QString::number( this->current_width_ * scale_factor ) );
	this->private_->ui_.y_axis_label_new_->setText( QString::fromUtf8("Y: ") + 
		QString::number( this->current_height_ * scale_factor ) );
	this->private_->ui_.z_axis_label_new_->setText( QString::fromUtf8("Z: ") + 
		QString::number( this->current_depth_ * scale_factor ) );
}
	
	
void LayerGroupWidget::enable_delete_button( bool enable )
{
	this->private_->ui_.delete_button_->setEnabled( enable );
}

void LayerGroupWidget::uncheck_delete_confirm()
{
    this->private_->ui_.confirm_delete_checkbox_->setChecked( false );
    this->private_->ui_.delete_button_->setEnabled( false );
    show_delete( false );
}

void LayerGroupWidget::show_layers( bool show )
{
	if( show ) 
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
	if( show ) 
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
		this->private_->ui_.group_resample_button_->setChecked( false );
	}
	show_selection_checkboxes( show );
	this->repaint();
}

void LayerGroupWidget::show_transform( bool show )
{
	if(show) 
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
		this->private_->ui_.group_transform_button_->setChecked( false );
	}
	show_selection_checkboxes( show );
	this->repaint();
}

void LayerGroupWidget::show_crop( bool show )
{
	if(show) 
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
		this->private_->ui_.group_crop_button_->setChecked( false );
	}
	show_selection_checkboxes( show );
	this->repaint();
}

void LayerGroupWidget::show_flip_rotate( bool show )
{
	if(show) 
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
		this->private_->ui_.group_flip_rotate_button_->setChecked( false );
		this->private_->ui_.flip_rotate_->hide();
	}
	show_selection_checkboxes( show );
	this->repaint();
}

void LayerGroupWidget::show_delete( bool show )
{
	if(show) 
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
		this->private_->ui_.group_delete_button_->setChecked( false );
	}
	show_selection_checkboxes( show );
	this->repaint();
}

	
}  //end namespace Seg3D






