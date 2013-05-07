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


#include <Application/BackscatterReconstruction/ViewerTool.h>
#include <Application/Tool/ToolFactory.h>

#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/Layer/LayerManager.h>

#include <Core/Viewer/Mouse.h>

#include <Core/Volume/DataVolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

// test
#include <iostream>
// test

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ViewerTool )

namespace Seg3D
{

ViewerTool::ViewerTool( const std::string& toolid ) :
SingleTargetTool( Core::VolumeType::MASK_E, toolid ) /*, private_( new ViewerToolPrivate )*/
{
}

ViewerTool::~ViewerTool()
{
  //this->disconnect_all();
}

void ViewerTool::execute( Core::ActionContextHandle context )
{
  // do nothing
}

void ViewerTool::activate()
{
  Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(),
                            ViewerManager::Instance()->layout_state_, ViewerManager::VIEW_1AND3_C );
  ViewerHandle viewer = ViewerManager::Instance()->get_viewer(0);
  viewer->view_mode_state_->set(Viewer::VOLUME_C);
}
  

}