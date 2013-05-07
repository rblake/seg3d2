/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#ifndef APPLICATION_BACKSCATTERRECONSTRUCTION_VIEWERTOOL_H
#define APPLICATION_BACKSCATTERRECONSTRUCTION_VIEWERTOOL_H

#include <Application/Tool/SingleTargetTool.h>
#include <Application/ViewerManager/ViewerManager.h>

namespace Seg3D
{
  
//class ViewerToolPrivate;
//typedef boost::shared_ptr< ViewerToolPrivate > ViewerToolPrivateHandle;

class ViewerTool : public SingleTargetTool
{
  
  SEG3D_TOOL
  (
   SEG3D_TOOL_NAME( "ViewerTool", "" )
   SEG3D_TOOL_MENULABEL( "Viewer" )
   SEG3D_TOOL_MENU( "Tools" )
   SEG3D_TOOL_SHORTCUT_KEY( "" )
   SEG3D_TOOL_URL( "" )
   SEG3D_TOOL_HOTKEYS( "" )
   SEG3D_TOOL_VERSION( "1.0" )
   )
  
public:
  ViewerTool( const std::string& toolid );
  virtual ~ViewerTool();
  
  virtual void execute( Core::ActionContextHandle context );    
  virtual void activate();    
};
  
}

#endif