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

#ifndef INTERFACE_APPINTERFACE_STYLESHEET_H
#define INTERFACE_APPINTERFACE_STYLESHEET_H

//Qt Includes
#include <QtGui>

namespace Seg3D
{
   
class StyleSheet
{

private:
	StyleSheet();
	
public:

	// === ViewerWidget Styles === //
		const static QString VIEWERWIDGET_TOOLBAR_C;
	
	// === StatusBar === //
		const static QString STATUSBAR_C;
		const static QString STATUSBAR_ERROR_C;
		const static QString STATUSBAR_WARNING_C;
		const static QString STATUSBAR_MESSAGE_C;
		const static QString STATUSBAR_DEBUG_C;
		
	// === ToolBoxPageWidget Styles === //
		// - page_background_ - //
		const static QString TOOLBOXPAGEWIDGET_PAGE_BACKGROUND_ACTIVE_C;
		const static QString TOOLBOXPAGEWIDGET_PAGE_BACKGROUND_INACTIVE_C;
		// - activate_button_ - //
		const static QString TOOLBOXPAGEWIDGET_ACTIVATE_BUTTON_ACTIVE_C;
		const static QString TOOLBOXPAGEWIDGET_ACTIVATE_BUTTON_INACTIVE_C;
		
	// === LayerImporterWidget === //
		// - Buttons - //
		const static QString LAYERIMPORTERWIDGET_ACTIVE_BUTTON_C;
		const static QString LAYERIMPORTERWIDGET_INACTIVE_BUTTON_C;

	
	// === LayerWidget Styles === //
		// - LayerWidget - //
		const static QString LAYERWIDGET_C;
		// - Volume Type Colors - //
		const static QString DATA_VOLUME_COLOR_C;
		const static QString MASK_VOLUME_COLOR_C;
		const static QString LABEL_VOLUME_COLOR_C;
		// - ACTIVE - //
		const static QString LAYER_WIDGET_BASE_ACTIVE_C;
		// - INACTIVE - //
		const static QString LAYER_WIDGET_BASE_INACTIVE_C;
		const static QString LAYER_WIDGET_LABEL_C;
		// - DROP - //
		const static QString LAYER_WIDGET_BASE_DROP_C;
		// - PICKED-UP - //
		const static QString LAYER_WIDGET_BASE_PICKED_UP_C;
		// - LOCKED - //
		const static QString LAYER_WIDGET_BASE_LOCKED_C;
		const static QString LAYER_WIDGET_BACKGROUND_LOCKED_C;
		const static QString LAYER_WIDGET_LABEL_LOCKED_C;
	
	// === LayerGroupWidget Styles === //
		// - ACTIVE - //
		const static QString GROUP_WIDGET_BASE_ACTIVE_C;
		const static QString GROUP_WIDGET_BACKGROUND_ACTIVE_C;
		const static QString GROUP_WIDGET_ACTIVATE_BUTTON_ACTIVE_C;
		// - INACTIVE - //
		const static QString GROUP_WIDGET_BASE_INACTIVE_C;
		const static QString GROUP_WIDGET_BACKGROUND_INACTIVE_C;
		const static QString GROUP_WIDGET_ACTIVATE_BUTTON_INACTIVE_C;
		
	// === PushDragButton's === //
		const static QString LAYER_PUSHDRAGBUTTON_C;
		const static QString GROUP_PUSHDRAGBUTTON_C;

	// === DropSpaceWidget === //
		const static QString DROPSPACEWIDGET_C;
	
	// === ViewerWidget === //
		const static QString VIEWERSELECTED_C;
		const static QString VIEWERNOTSELECTED_C;

};    

}  // end namespace Seg3D

#endif