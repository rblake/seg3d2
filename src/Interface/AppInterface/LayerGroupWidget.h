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

#ifndef INTERFACE_APPINTERFACE_LAYERGROUPWIDGET_H
#define INTERFACE_APPINTERFACE_LAYERGROUPWIDGET_H

// QT Includes
#include <QtGui>

// Application includes
#include <Application/Layer/LayerGroup.h>
#include <Interface/AppInterface/LayerWidget.h>

namespace Seg3D
{

class LayerGroupWidgetPrivate;

class LayerGroupWidget : public QWidget
{
	Q_OBJECT
	
Q_SIGNALS:
	// PREP_GROUPS_FOR_DRAG_AND_DROP:
	// this signal tells the LayerManagerWidget that it needs to prep all the groups for drag and 
	// drop.
	void prep_layers_for_drag_and_drop_signal_( bool );
	
	void prep_groups_for_drag_and_drop_signal_( bool );
	
	void picked_up_group_size_signal_( int );
	
	// -- constructor/destructor --
public:
	LayerGroupWidget( QWidget* parent, LayerHandle layer );
	virtual ~LayerGroupWidget();
	
public Q_SLOTS:
	void show_layers( bool show );
	void show_resample( bool show );
	void show_transform( bool show );
	void show_crop( bool show );
	void show_flip_rotate( bool show );
	void show_delete( bool show );
	void show_selection_checkboxes( bool show );
	void enable_delete_button( bool enable );
	
	// NOTIFY_LAYER_MANAGER_WIDGET:
	// this function triggers the prep_groups_for_drag_and_drop signal.
	void notify_layer_manager_widget( bool move_time );
	
	// PREP_FOR_ANIMATION:
	// this function replaces the widget with a screenshot of the widget for speed
	void prep_for_animation( bool move_time );
	
	
public:
	//void add_layer( LayerHandle layer );
	void insert_layer( LayerHandle layer, int index );
	bool delete_layer( LayerHandle layer );
	const std::string &get_group_id();
	void set_active( bool active );
	LayerWidgetQWeakHandle set_active_layer( LayerHandle layer );
	void seethrough( bool see );
	
	// PREP_LAYERS_FOR_DRAG_AND_DROP:
	// this function tells each layer to pepare for drag and drop by replacing the actual widgets
	// with images of themselves.
	void prep_layers_for_drag_and_drop( bool move_time );
	
	void set_picked_up_group_size( int group_height );


	
	
protected:
	void resizeEvent( QResizeEvent *event );
	
private Q_SLOTS:
    void adjust_new_size_labels( double scale_factor );
    void uncheck_delete_confirm();
    void set_picked_up( bool up ){ this->picked_up_ = up; }
    void set_drop( bool drop );
   
	void set_drop_target( LayerGroupWidget* target_layer );
	void mousePressEvent( QMouseEvent* event );
	void dropEvent( QDropEvent* event );
	void dragEnterEvent( QDragEnterEvent* event );
	void dragLeaveEvent( QDragLeaveEvent* event );
	void show_group();
	void hide_group();

	// -- widget internals --
private:
    boost::shared_ptr< LayerGroupWidgetPrivate > private_;
    QVector< LayerWidgetQHandle > layer_list_;
    
private:
	std::string group_id_;
	int picked_up_group_height_;
	int current_height_;
	int current_width_;
	int current_depth_;
	bool group_menus_open_;
	bool picked_up_;
	bool drop_group_set_;
	LayerGroupWidget* drop_group_;
		
};
		


} //end namespace Seg3D

#endif