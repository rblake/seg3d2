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
	
public:
	void add_layer( LayerHandle layer );
	void delete_layer( LayerHandle layer );
	std::string &get_group_id();
	void set_active( bool active );
	void set_active_layer( LayerHandle layer );
	void clear_all_layers()
    {
        this->layer_list_.clear();
    }

private Q_SLOTS:
    void adjust_new_size_labels( double scale_factor );
    void uncheck_delete_confirm();
    
private:

	// Drag and drop functions
	void dropEvent(QDropEvent* event);
	void mousePressEvent(QMouseEvent *event);
	void dragEnterEvent(QDragEnterEvent* event);
	void dragLeaveEvent( QDragLeaveEvent * event );
	
	// Validate that the event location is a LayerWidget and return a handle to it
	LayerWidget_handle validate_location( const QPoint& point );
	
	// Check to see if the layer_id is currently in the group
	LayerWidget_handle check_for_layer( const std::string &layer );

	
	// -- widget internals --
private:
    boost::shared_ptr< LayerGroupWidgetPrivate > private_;
    QVector< LayerWidget_handle > layer_list_;
			
};
		


} //end namespace Seg3D

#endif