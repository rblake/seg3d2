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

#ifndef APPLICATION_BACKSCATTERRECONSTRUCTION_CALIBRATIONTOOL_H
#define APPLICATION_BACKSCATTERRECONSTRUCTION_CALIBRATIONTOOL_H

#include <Application/Tool/SingleTargetTool.h>

namespace Seg3D
{

class CalibrationToolPrivate;
typedef boost::shared_ptr< CalibrationToolPrivate > CalibrationToolPrivateHandle;

class CalibrationTool : public SingleTargetTool
{

SEG3D_TOOL
(
  SEG3D_TOOL_NAME( "CalibrationTool", "" )
  SEG3D_TOOL_MENULABEL( "Calibration" )
  SEG3D_TOOL_MENU( "Tools" )
  SEG3D_TOOL_SHORTCUT_KEY( "" )
  SEG3D_TOOL_URL( "" )
  SEG3D_TOOL_HOTKEYS( "" )
  SEG3D_TOOL_VERSION( "1.0" )
)

public:
  CalibrationTool( const std::string& toolid );
  virtual ~CalibrationTool();

	virtual void execute( Core::ActionContextHandle context );

	/// ACTIVATE:
	/// Activate a tool: this tool is set as the active tool and hence it should
	/// setup the right mouse tools in the viewers.
	virtual void activate();

  void save( Core::ActionContextHandle context, int index, std::string layerid );
	void segment( Core::ActionContextHandle context );
  
	// -- state --
  // target_layer_state_ is selected layer

  Core::StateStringHandle input_data_id_;
  Core::StateStringVectorHandle calibrationSet_state_;
  Core::StateStringVectorHandle maskLayers_state_;

private:
	CalibrationToolPrivateHandle private_;
};
  
}

#endif