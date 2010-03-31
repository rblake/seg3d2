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
#include <QtGui>

// Application Includes
#include <Utils/Core/EnumClass.h>

#include <Application/Layer/LayerGroup.h>
#include <Interface/AppInterface/LayerGroupWidget.h>

namespace Seg3D
{
	
// enum for layer types
SCI_ENUM_CLASS
(
	LayerType,
	DATA_LAYER_E = 1, 
	MASK_LAYER_E = 2, 
	LABEL_LAYER_E = 3
)
	
	

//class LayerManagerWidgetPrivate;
//typedef boost::shared_ptr< LayerManagerWidgetPrivate > LayerManagerWidgetPrivateHandle_type;

class LayerManagerWidget : public QScrollArea
{
	// Needed to make it a Qt object
Q_OBJECT

//constructor - destructor
public:
	LayerManagerWidget( QWidget *parent = 0 );
	virtual ~LayerManagerWidget();


public:
	void insert_layer( LayerHandle layer );
	void delete_layer( LayerGroupHandle group );
	bool refresh_group( LayerGroupHandle group );
	void make_new_group( LayerHandle layer );
	void delete_group( LayerGroupHandle group );
	void show_group( LayerGroupHandle group );
	void set_active_group( LayerGroupHandle group );
	void set_active_layer( LayerHandle layer );
	

private:
	// private Qt GUI Components for the LayerManagerWidget
	QWidget* main_;
	QVBoxLayout* main_layout_;
	QVBoxLayout* group_layout_;
	


private:
	QList< LayerGroupWidget_handle > group_list_;
};

} //endnamespace Seg3d

#endif
