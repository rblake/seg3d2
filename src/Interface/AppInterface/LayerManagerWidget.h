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

#ifndef INTERFACE_APPINTERFACE_LAYERMANAGERWIDGET_H
#define INTERFACE_APPINTERFACE_LAYERMANAGERWIDGET_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// Qt includes
#include <QtGui/QWidget>
#include <QtGui/QScrollArea>

// Core includes
#include <Core/Utils/EnumClass.h>
#include <Core/Utils/ConnectionHandler.h>

// Application includes
#include <Application/Layer/LayerGroup.h>

// Interface includes
#include <Interface/AppInterface/LayerGroupWidget.h>

namespace Seg3D
{


// Forward declarations	
class LayerManagerWidgetPrivate;
typedef boost::shared_ptr<LayerManagerWidgetPrivate> LayerManagerWidgetPrivateHandle;
	
class LayerManagerWidget : public QScrollArea, public Core::ConnectionHandler
{
// Need to make it a Qt object
Q_OBJECT
	
	// -- constructor/destructor
public:
	LayerManagerWidget( QWidget *parent = 0 );
	virtual ~LayerManagerWidget();

	// -- function that handle signals --
private:
	// HANDLE_GROUP_INTERNALS_CHANGE:
	// This function is called when the internals of a group are modified.
	void handle_group_internals_change( LayerGroupHandle group );
	
	// HANDLE_GROUP_CHANGED:
	// This function is called when new groups are added or deleted and when their
	// order changes.
	void handle_groups_changed();

	// PRE_LOAD_STATES:
	// Called on pre_load_states_signal_ from the StateEngine.
	void pre_load_states();

	// POST_LOAD_STATES:
	// Called on post_load_states_signal_ from the StateEngine.
	void post_load_states();

	// RESET:
	// Remove all the layer widgets.
	void reset();

	// MAKE_NEW_GROUP:
	// function that creates a new group to put layers into. 
	LayerGroupWidget* make_new_group( LayerGroupHandle group );
	
	
private Q_SLOTS:
	// PREP_LAYERS_FOR_DRAG_AND_DROP:
	// this function tells the groups to prepare their layers for drag and drop
	void prep_layers_for_drag_and_drop( bool move_time );
	
	// PREP_GROUPS_FOR_DRAG_AND_DROP:
	// this function tells the groups to prepare for drag and drop
	void prep_groups_for_drag_and_drop( bool move_time );
	
	// NOTIFY_PICKED_UP_GROUP_SIZE:
	// this function will notify all the groups of what size the currently picked up group is
	void notify_picked_up_group_size( int group_size );
	
	void notify_groups_of_picked_up_layer_size( int layer_size );
	
	// -- internals of this class --
private:
	friend class LayerManagerWidgetPrivate;

	LayerManagerWidgetPrivateHandle private_;
};

} //endnamespace Seg3d

#endif
