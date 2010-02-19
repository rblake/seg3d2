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

#include <Interface/AppInterface/LayerManagerWidget.h>
#include <Utils/Core/Log.h>

#include <sstream>
#include <iostream>
#include <boost/lexical_cast.hpp>

//Interface Includes
#include <Interface/ToolInterface/CustomWidgets/SliderSpinComboInt.h>
#include <Interface/ToolInterface/CustomWidgets/SliderSpinComboDouble.h>
#include <Interface/QtInterface/QtBridge.h>

// Gt Gui Includes
#include "ui_LayerGroupWidget.h"
#include "ui_LayerWidget.h"

namespace Seg3D  {
	
  typedef boost::shared_ptr <QWidget> QWidgetHandle;
  typedef std::list<QWidgetHandle> newlayerlisttype;

typedef boost::shared_ptr <Ui::LayerGroupWidget> LayerGroupWidgetHandle_type;
typedef std::list<LayerGroupWidgetHandle_type> GroupList_type;

typedef boost::shared_ptr <Ui::LayerWidget> LayerWidgetHandle_type;
typedef std::list<LayerWidgetHandle_type> LayerList_type;


class LayerManagerWidgetPrivate {
public:
	Ui::LayerGroupWidget group_ui_;
	Ui::LayerWidget layer_ui_;
	
	GroupList_type group_list_;
	LayerList_type layer_list_;
  newlayerlisttype layer_widget_list_;
};


LayerManagerWidget::LayerManagerWidget( QWidget* parent ) :
QScrollArea( parent )
{ 
	//initialize the number of groups
	number_of_groups_ = 0;
	
	// make a new LayerManagerWidgetPrivateHandle_type object
	private_ = LayerManagerWidgetPrivateHandle_type ( new LayerManagerWidgetPrivate );
	
	// set some values for the scrollarea widget
	setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
	setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded );
	setContentsMargins( 1, 1, 1, 1 );
	setWidgetResizable( true );  
	
	main_ = new QWidget( parent );
	setWidget( main_ );
	
	main_layout_ = new QVBoxLayout( main_ );
	main_layout_->setContentsMargins( 1, 1, 1, 1 );
	main_layout_->setSpacing( 1 );
	
	group_layout_ = new QVBoxLayout;
	group_layout_->setSpacing(5);
	main_layout_->addLayout( group_layout_ );
	main_layout_->addStretch();
	
	main_->setLayout( main_layout_ );
	main_->setSizePolicy( QSizePolicy::MinimumExpanding, QSizePolicy::Preferred );
	
	
	// add test groups
	new_group( QString::fromUtf8( "600x300x240" ));
	new_group( QString::fromUtf8( "600x300x240" ));
	new_group( QString::fromUtf8( "400x300x340" ));
	new_group( QString::fromUtf8( "300x200x400" ));
	new_group( QString::fromUtf8( "800x300x940" ));
	new_group( QString::fromUtf8( "900x200x300" ));
	new_group( QString::fromUtf8( "700x300x440" ));
	new_group( QString::fromUtf8( "100x200x200" ));
	
	//add test layers
	new_layer(DATA_LAYER_E, QString::fromUtf8( "Mouse Study 23430.43" ), QString::fromUtf8( "300x200x400" ));
	new_layer(DATA_LAYER_E, QString::fromUtf8( "Rat Spine-432.d9s" ), QString::fromUtf8("300x200x400") );
	new_layer(DATA_LAYER_E, QString::fromUtf8( "CervicalSpine_23k3s.23" ), QString::fromUtf8("300x200x400") );
	new_layer(DATA_LAYER_E, QString::fromUtf8( "Neuro Study Rat.3421" ), QString::fromUtf8("300x200x400") );
	new_layer(DATA_LAYER_E, QString::fromUtf8( "Cervical_239.sle" ), QString::fromUtf8("700x300x440") );
	
}
// destructor
LayerManagerWidget::~LayerManagerWidget()
{ }
	
	
	// create a new group and add it to the group_layout_
	// add it also to the group_list_
	void LayerManagerWidget::new_group( const QString &dimensions )
	{
		if ( !validate_new_group(dimensions) )
		{
			return;
		}
		
		// create a new widget to add to the group_layout_
		QWidget* new_group_ = new QWidget();
		
		private_->group_ui_.setupUi( new_group_ );
		
		// create a handle to the group ui object
		LayerGroupWidgetHandle_type group_handle_ ( new Ui::LayerGroupWidget( private_->group_ui_ ));
		
		// set some values of the GUI
		group_handle_->activate_button_->setText( dimensions );
		group_handle_->group_frame_->hide();
		group_handle_->group_tools_->hide();
		
		// hide the tool bars
		group_handle_->roi_->hide();
		group_handle_->resample_->hide();
		group_handle_->flip_rotate_->hide();
		group_handle_->transform_->hide();
		
		// add the slider spinner combo's for the crop
		SliderSpinComboInt* x_adjuster_crop = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_11->addWidget(x_adjuster_crop);
		x_adjuster_crop->setObjectName(QString::fromUtf8("x_adjuster_crop"));
		
		SliderSpinComboInt* y_adjuster_crop = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_12->addWidget(y_adjuster_crop);
		y_adjuster_crop->setObjectName(QString::fromUtf8("y_adjuster_crop"));
		
		SliderSpinComboInt* z_adjuster_crop = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_14->addWidget(z_adjuster_crop);
		z_adjuster_crop->setObjectName(QString::fromUtf8("z_adjuster_crop"));
		
		SliderSpinComboInt* height_adjuster_crop = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_7->addWidget(height_adjuster_crop);
		height_adjuster_crop->setObjectName(QString::fromUtf8("height_adjuster_crop"));
		
		SliderSpinComboInt* width_adjuster_crop = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_9->addWidget(width_adjuster_crop);
		width_adjuster_crop->setObjectName(QString::fromUtf8("width_adjuster_crop"));
		
		SliderSpinComboInt* depth_adjuster_crop = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_10->addWidget(depth_adjuster_crop);
		depth_adjuster_crop->setObjectName(QString::fromUtf8("depth_adjuster_crop"));
		
		
		// add the slider spinner combo's for the tranform
		SliderSpinComboInt* x_adjuster_transform = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_16->addWidget(x_adjuster_transform);
		x_adjuster_transform->setObjectName(QString::fromUtf8("x_adjuster_transform"));
		
		SliderSpinComboInt* y_adjuster_transform = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_17->addWidget(y_adjuster_transform);
		y_adjuster_transform->setObjectName(QString::fromUtf8("y_adjuster_transform"));
		
		SliderSpinComboInt* z_adjuster_transform = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_18->addWidget(z_adjuster_transform);
		z_adjuster_transform->setObjectName(QString::fromUtf8("z_adjuster_transform"));
		
		SliderSpinComboDouble* height_adjuster_transform = new SliderSpinComboDouble(new_group_);
		group_handle_->horizontalLayout_20->addWidget(height_adjuster_transform);
		height_adjuster_transform->setObjectName(QString::fromUtf8("height_adjuster_transform"));
		
		SliderSpinComboDouble* width_adjuster_transform = new SliderSpinComboDouble(new_group_);
		group_handle_->horizontalLayout_19->addWidget(width_adjuster_transform);
		width_adjuster_transform->setObjectName(QString::fromUtf8("width_adjuster_transform"));
		
		SliderSpinComboDouble* depth_adjuster_transform = new SliderSpinComboDouble(new_group_);
		group_handle_->horizontalLayout_21->addWidget(depth_adjuster_transform);
		depth_adjuster_transform->setObjectName(QString::fromUtf8("depth_adjuster_transform"));
		
		
		SliderSpinComboInt* scale_adjuster = new SliderSpinComboInt(new_group_);
		group_handle_->horizontalLayout_15->addWidget(scale_adjuster);
		scale_adjuster->setObjectName(QString::fromUtf8("scale_adjuster"));
		
		
		// connect the signals and slots
		connect(group_handle_->open_button_, SIGNAL( clicked( bool )), this, SLOT(hide_show_layers( bool )));
		connect(group_handle_->group_resample_button_, SIGNAL( clicked( bool ) ), this, SLOT(hide_show_resample( bool )));
		connect(group_handle_->group_crop_button_, SIGNAL( clicked( bool )), this, SLOT(hide_show_roi( bool )));
		connect(group_handle_->group_transform_button_, SIGNAL( clicked( bool ) ), this, SLOT(hide_show_transform( bool )));
		connect(group_handle_->group_flip_rotate_button_, SIGNAL( clicked ( bool ) ), this, SLOT(hide_show_flip_rotate( bool )));
		
		// add the new widget to the group_layout_
		group_layout_->addWidget( new_group_ );
		
		// add the new group handle to the group_list_
		private_->group_list_.push_back(group_handle_);
		
		// increment the number of groups
		number_of_groups_++;
	} // end new_group
	
	
	// iterate through all the current layers and check to see if a layer with similar dimensions already exists
	bool LayerManagerWidget::validate_new_group( const QString &dimensions )
	{
		for ( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++ )
		{
			if ( (*i)->activate_button_->text() == dimensions )
			{ 
				SCI_LOG_ERROR(std::string( "A group with dimensions: '") +  dimensions.toStdString() + "' already exists, not creating another." );
				return false; 
			}
		}
		return true;
	}
	
	
	// create a new layer and add it to the proper group
	// add it also to the group_list_
	void LayerManagerWidget::new_layer( int type, const QString &label, const QString &dimensions )
	{

		// create a new widget to add to the group layout
		QWidget* new_layer_ = new QWidget();
		private_->layer_ui_.setupUi( new_layer_ );
		
		// create a handle to the layer ui object
		LayerWidgetHandle_type layer_handle_ ( new Ui::LayerWidget( private_->layer_ui_ ));
    
    //QWidgetHandle layer_widget_handle_ (new QWidget(new_layer_));

		
		// set some values in the GUI
		layer_handle_->label_->setText(label);
    //->findChild< QToolButton* >( "close_button_" )
    //SCI_LOG_ERROR(new_layer_->findChild< QLineEdit* >( "label_" )->text().toStdString());
    //SCI_LOG_ERROR(layer_widget_handle_->findChild< QLineEdit* >( "label_" )->text().toStdString());//->setText(QString::fromUtf8("testing"));
    
		layer_handle_->dimensions_->setText(dimensions);
		
		// hide the tool bars and the selection checkbox
		layer_handle_->color_bar_->hide();
		layer_handle_->bright_contrast_bar_->hide();
		layer_handle_->checkbox_widget_->hide();
		layer_handle_->opacity_bar_->hide();
		layer_handle_->progress_bar_bar_->hide();
    layer_handle_->border_bar_->hide();
		layer_handle_->dimensions_->hide();
		
		SliderSpinComboInt* opacity_adjuster_ = new SliderSpinComboInt(new_layer_);
		layer_handle_->verticalLayout_2->addWidget(opacity_adjuster_);
		opacity_adjuster_->setObjectName(QString::fromUtf8("opacity_adjuster_"));
		
		SliderSpinComboInt* brightness_adjuster_ = new SliderSpinComboInt(new_layer_);
		layer_handle_->brightness_h_layout_->addWidget(brightness_adjuster_);
		brightness_adjuster_->setObjectName(QString::fromUtf8("brightness_adjuster_"));
		
		SliderSpinComboInt* contrast_adjuster_ = new SliderSpinComboInt(new_layer_);
		layer_handle_->contrast_h_layout_->addWidget(contrast_adjuster_);
		contrast_adjuster_->setObjectName(QString::fromUtf8("contrast_adjuster_"));

    layer_handle_->border_selection_combo_->addItem(QString::fromUtf8("No Border"));
    layer_handle_->border_selection_combo_->addItem(QString::fromUtf8("Thin Border"));
    layer_handle_->border_selection_combo_->addItem(QString::fromUtf8("Thick Border"));
		
		
		// connect the signals and slots
		connect( layer_handle_->color_button_,					SIGNAL( clicked( bool )),	this, SLOT( hide_show_color_choose_bar( bool )));
		connect( layer_handle_->opacity_button_,				SIGNAL( clicked( bool )),	this, SLOT( hide_show_opacity_bar( bool )));
		connect( layer_handle_->brightness_contrast_button_,	SIGNAL( clicked( bool )),	this, SLOT( hide_show_brightness_contrast_bar( bool )));
		connect( layer_handle_->compute_iso_surface_button_,	SIGNAL( clicked() ),			this, SLOT( show_progress_bar_bar()));
		connect( layer_handle_->progress_cancel_button_,		SIGNAL( clicked() ),			this, SLOT( hide_progress_bar_bar()));
		connect( layer_handle_->lock_button_,					SIGNAL( clicked( bool )),	this, SLOT( lock_unlock_layer( bool )));
    connect( layer_handle_->fill_border_button_,  SIGNAL( clicked( bool )), this, SLOT( hide_show_border_bar ( bool )));
		
		// add the new widget to the appropriate group's group_frame_layout_
		for ( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++ )
		{
			if ( (*i)->activate_button_->text() == dimensions )
			{
				(*i)->group_frame_layout_->addWidget( new_layer_ );
				(*i)->group_frame_->show();
				(*i)->group_tools_->show();
				(*i)->open_button_->setChecked( true );
				(*i)->group_visibility_button_->setChecked( true );
				break;
			}
		}
		
		// add the new layer handle to the layer_list_
		private_->layer_list_.push_back(layer_handle_);
		
	} // end new_layer
    
	
	// --- Functions for Hiding and showing Layer and Group options --- //
	
	// hide or show the Resample tool box
	void LayerManagerWidget::hide_show_resample( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		QString group_dimensions_;
		
		for( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++ )
		{
			if( (*i)->group_resample_button_ == hide_show_button )
			{
				group_dimensions_ = (*i)->activate_button_->text();
				if( hideshow ) 
				{
					(*i)->resample_->show();
					(*i)->roi_->hide();
					(*i)->group_crop_button_->setChecked(false);
					(*i)->flip_rotate_->hide();
					(*i)->group_flip_rotate_button_->setChecked( false );
					(*i)->transform_->hide();
					(*i)->group_transform_button_->setChecked( false );
				}
				else 
				{
					(*i)->resample_->hide();
				}
			} 
		}
		
		// hide or show the appropriate checkboxes
		turn_off_or_on_checkboxes(group_dimensions_, hideshow);
		
	} // end hide_show_resample
	
	
	// hide or show the Region Of Interest tool box
	void LayerManagerWidget::hide_show_roi( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		QString group_dimensions_;
		
		for ( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++ )
		{
			if ( (*i)->group_crop_button_ == hide_show_button )
			{
				group_dimensions_ = (*i)->activate_button_->text();
				if ( hideshow ) 
				{
					(*i)->roi_->show();
					(*i)->resample_->hide();
					(*i)->group_resample_button_->setChecked( false );
					(*i)->flip_rotate_->hide();
					(*i)->group_flip_rotate_button_->setChecked( false );
					(*i)->transform_->hide();
					(*i)->group_transform_button_->setChecked( false );
				}
				else
				{
					(*i)->roi_->hide();
				}
			} 
		} 
		
		// hide or show the appropriate checkboxes
		turn_off_or_on_checkboxes(group_dimensions_, hideshow);
		
	} // end hide_show_roi
	
	
	// hide or show the Region Of Interest tool box
	void LayerManagerWidget::hide_show_transform( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		QString group_dimensions_;
		
		for ( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++ )
		{
			if ( (*i)->group_transform_button_ == hide_show_button )
			{
				group_dimensions_ = (*i)->activate_button_->text();
				if ( hideshow ) 
				{
					(*i)->transform_->show();
					(*i)->resample_->hide();
					(*i)->roi_->hide();
					(*i)->flip_rotate_->hide();
					(*i)->group_flip_rotate_button_->setChecked( false );
					(*i)->group_resample_button_->setChecked( false );
					(*i)->group_crop_button_->setChecked( false );
				}
				else
				{
					(*i)->transform_->hide();
				}
			} 
		} 
		// hide or show the appropriate checkboxes
		turn_off_or_on_checkboxes(group_dimensions_, hideshow);

	} // end hide_show_transform
	
	
	// hide or show the flip rotate tool box
	void LayerManagerWidget::hide_show_flip_rotate( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		QString group_dimensions_;
		
		for ( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++ )
		{
			if ( (*i)->group_flip_rotate_button_ == hide_show_button )
			{
				group_dimensions_ = (*i)->activate_button_->text();
				if ( hideshow ) 
				{
					(*i)->flip_rotate_->show();
					(*i)->transform_->hide();
					(*i)->resample_->hide();
					(*i)->roi_->hide();
					(*i)->group_resample_button_->setChecked( false );
					(*i)->group_crop_button_->setChecked( false );
					(*i)->group_transform_button_->setChecked( false );
				}
				else
				{
					(*i)->flip_rotate_->hide();
				}
			} 
		} 
		// hide or show the appropriate checkboxes
		turn_off_or_on_checkboxes(group_dimensions_, hideshow);
		
	} // end hide_show_flip_rotate
	
	
	// turn off or on the checkboxes
	void LayerManagerWidget::turn_off_or_on_checkboxes( const QString &dimensions, const bool &hideshow )
	{
		for ( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( hideshow && ( dimensions == (*i)->dimensions_->text())) 
			{
				if( !(*i)->lock_button_->isChecked() )
					(*i)->checkbox_widget_->show();
			}
			else {
				if( dimensions == (*i)->dimensions_->text())
				{
					(*i)->checkbox_widget_->hide();
				}
			}
		}
	}
	
	
	
	// hide or show the layers
	void LayerManagerWidget::hide_show_layers( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( GroupList_type::const_iterator i = private_->group_list_.begin(); i != private_->group_list_.end(); i++)
		{
			if( (*i)->open_button_ == hide_show_button )
			{
				if( hideshow ) 
				{
					(*i)->group_frame_->show();
					(*i)->group_tools_->show();
				}
				else 
				{
					(*i)->group_frame_->hide();
					(*i)->group_tools_->hide();
				}
			} 
		}
	} // end hide_show_layers
	
	
	
	// hide or show the color choose bar
	void LayerManagerWidget::hide_show_color_choose_bar( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->color_button_ == hide_show_button )
			{
				if( hideshow ) 
				{
					(*i)->color_bar_->show();
					(*i)->bright_contrast_bar_->hide();
					(*i)->brightness_contrast_button_->setChecked( false );
					(*i)->opacity_bar_->hide();
					(*i)->opacity_button_->setChecked( false );
          (*i)->border_bar_->hide();
					(*i)->fill_border_button_->setChecked( false );
					//TEST CODE
					(*i)->progress_bar_bar_->hide();
					(*i)->compute_iso_surface_button_->setChecked( false );
				}
				else
				{
					(*i)->color_bar_->hide();
				}
			} 
		}
	} // end hide_show_color_choose_bar
	
	
	
	
	// hide or show the opacity bar
	void LayerManagerWidget::hide_show_opacity_bar( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->opacity_button_ == hide_show_button )
			{
				if( hideshow ) 
				{
					(*i)->opacity_bar_->show();
					(*i)->bright_contrast_bar_->hide();
					(*i)->brightness_contrast_button_->setChecked( false );
					(*i)->color_bar_->hide();
					(*i)->color_button_->setChecked( false );
          (*i)->border_bar_->hide();
					(*i)->fill_border_button_->setChecked( false );
					//TEST CODE
					(*i)->progress_bar_bar_->hide();
					(*i)->compute_iso_surface_button_->setChecked( false );
				}
				else
				{
					(*i)->opacity_bar_->hide();
				}
			} 
		}
	} // end hide_show_opacity_bar
	
	
	// hide or show the brightness/contrast bar
	void LayerManagerWidget::hide_show_brightness_contrast_bar( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->brightness_contrast_button_ == hide_show_button )
			{
				if( hideshow ) 
				{
					(*i)->bright_contrast_bar_->show();
					(*i)->color_bar_->hide();
					(*i)->color_button_->setChecked( false );
					(*i)->opacity_bar_->hide();
					(*i)->opacity_button_->setChecked( false );
          (*i)->border_bar_->hide();
					(*i)->fill_border_button_->setChecked( false );
					//TEST CODE
					(*i)->progress_bar_bar_->hide();
					(*i)->compute_iso_surface_button_->setChecked( false );
				}
				else
				{
					(*i)->bright_contrast_bar_->hide();
				}
			} 
		}
	} // end hide_show_brightness_contrast_bar


  // hide or show the brightness/contrast bar
	void LayerManagerWidget::hide_show_border_bar( bool hideshow )
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->fill_border_button_ == hide_show_button )
			{
				if( hideshow ) 
				{
					(*i)->border_bar_->show();
          (*i)->bright_contrast_bar_->hide();
					(*i)->brightness_contrast_button_->setChecked( false );
					(*i)->color_bar_->hide();
					(*i)->color_button_->setChecked( false );
					(*i)->opacity_bar_->hide();
					(*i)->opacity_button_->setChecked( false );
					//TEST CODE
					(*i)->progress_bar_bar_->hide();
					(*i)->compute_iso_surface_button_->setChecked( false );
				}
				else
				{
					(*i)->border_bar_->hide();
				}
			} 
		}
	} // end hide_show_brightness_contrast_bar


	
	
	// show the progress bar bar
	void LayerManagerWidget::show_progress_bar_bar()
	{
		QToolButton *hide_show_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->compute_iso_surface_button_ == hide_show_button )
			{
				(*i)->progress_bar_bar_->show();
			}
		}
	} // end show_progress_bar_bar
	
	// hide the progress bar bar
	void LayerManagerWidget::hide_progress_bar_bar()
	{
		QToolButton *cancel_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->progress_cancel_button_ == cancel_button )
			{
				(*i)->progress_bar_bar_->hide();
			}
		}
	} // end hide_progress_bar_bar
	
	
	
	void LayerManagerWidget::color_button_clicked()
	{
		QToolButton *color_button = ::qobject_cast<QToolButton*>( sender() );
		
		for ( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if ( (*i)->color_button_01_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(251,255,74);}" ));
			}
			if ( (*i)->color_button_02_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(248,188,37);}" ));
			}
			if ( (*i)->color_button_03_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(249,152,28);}" ));
			}
			if ( (*i)->color_button_04_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(251,78,15);}" ));
			}
			if ( (*i)->color_button_05_ == color_button ) 
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(252,21,17);}" ));
			}
			if ( (*i)->color_button_06_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(166,12,73);}" ));
			}
			if ( (*i)->color_button_07_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(135,0,172);}" ));
			}
			if ( (*i)->color_button_08_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(66,0,161);}" ));
			}
			if ( (*i)->color_button_09_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(45,66,250);}" ));
			}
			if ( (*i)->color_button_10_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(32,146,204);}" ));
			}
			if ( (*i)->color_button_11_ == color_button ) 
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(100,177,61);}" ));
			}
			if ( (*i)->color_button_12_ == color_button )
			{
				(*i)->typeBackground_->setStyleSheet( QString::fromUtf8( "QWidget#typeBackground_{background-color: rgb(205,235,66);}" ));
			}
		}
	} // end color_button_clicked


  // lock or unlock a layer
	void LayerManagerWidget::lock_unlock_layer( bool lockunlock )
	{
		QToolButton *lock_unlock_button = ::qobject_cast<QToolButton*>( sender() );
		
		for( LayerList_type::const_iterator i = private_->layer_list_.begin(); i != private_->layer_list_.end(); i++ )
		{
			if( (*i)->lock_button_ == lock_unlock_button )
			{
				if( lockunlock ) 
				{
					(*i)->header_->setStyleSheet(QString::fromUtf8("QWidget#header_{"
					  "background-color: qlineargradient(spread:pad, x1:0.5, y1:0, x2:0.5, y2:0.960227, stop:0 rgba(201, 201, 201, 255), stop:0.155779 rgba(208, 208, 208, 255), stop:1 rgba(184, 184, 184, 255));"
					  "border-radius: 6px;}"));
					(*i)->typeBackground_->setStyleSheet(QString::fromUtf8(
							  "QWidget#typeBackground_{"
							  "background-color: gray;"
							  "border: 1px solid rgb(141, 141, 141);"
							  "border-radius: 4px;}"));
					
					(*i)->label_->setStyleSheet(QString::fromUtf8("QLineEdit#label_{background-color:rgb(208, 208, 208); color: gray;}"));
					(*i)->colorChooseButton_->setEnabled( false );
					(*i)->opacity_button_->setEnabled( false );
					(*i)->visibility_button_->setEnabled( false );
					(*i)->color_button_->setEnabled( false );
					(*i)->compute_iso_surface_button_->setEnabled( false );
					(*i)->iso_surface_button_->setEnabled( false );
					(*i)->fill_border_button_->setEnabled( false );
					(*i)->volume_rendered_button_->setEnabled( false );
					(*i)->brightness_contrast_button_->setEnabled( false );
					(*i)->checkbox_widget_->hide();

				}
				else
				{
					(*i)->header_->setStyleSheet(QString::fromUtf8(
						"QWidget#header_{"
						"background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0.512563 rgba(255, 255, 255, 0));}"));
					(*i)->typeBackground_->setStyleSheet(QString::fromUtf8(
						"QWidget#typeBackground_{"
					    "background-color: rgb(255, 128, 0);"
					    "border: 1px solid rgb(141, 141, 141);"
					    "border-radius: 4px;}"));
					(*i)->label_->setStyleSheet(QString::fromUtf8(
						"QLineEdit#label_{"
						"	text-align: left;"
						"	color: black;"
						"	margin-right: 3px;"
						"	background-color: rgb(243, 243, 243);}"));

					(*i)->colorChooseButton_->setEnabled( true );
					(*i)->opacity_button_->setEnabled( true );
					(*i)->visibility_button_->setEnabled( true );
					(*i)->color_button_->setEnabled( true );
					(*i)->compute_iso_surface_button_->setEnabled( true );
					(*i)->iso_surface_button_->setEnabled( true );
					(*i)->fill_border_button_->setEnabled( true );
					(*i)->volume_rendered_button_->setEnabled( true );
					(*i)->brightness_contrast_button_->setEnabled( true );
					
				}
			} 
		}
	} // end lock_unlock_layer
	
	
	
} //end Seg3D namespace
