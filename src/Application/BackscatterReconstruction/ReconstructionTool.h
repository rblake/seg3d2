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

#ifndef APPLICATION_BACKSCATTERRECONSTRUCTION_RECONSTRUCTIONTOOL_H
#define APPLICATION_BACKSCATTERRECONSTRUCTION_RECONSTRUCTIONTOOL_H

#include <Application/Tool/SingleTargetTool.h>

namespace Seg3D
{
  
class ReconstructionToolPrivate;
typedef boost::shared_ptr< ReconstructionToolPrivate > ReconstructionToolPrivateHandle;

class ReconstructionTool : public SingleTargetTool
{
  
SEG3D_TOOL
(
  SEG3D_TOOL_NAME( "ReconstructionTool", "" )
  SEG3D_TOOL_MENULABEL( "Reconstruction" )
  SEG3D_TOOL_MENU( "Tools" )
  SEG3D_TOOL_SHORTCUT_KEY( "" )
  SEG3D_TOOL_URL( "" )
  SEG3D_TOOL_HOTKEYS( "" )
  SEG3D_TOOL_VERSION( "1.0" )
)
  
public:
  ReconstructionTool( const std::string& toolid );
  virtual ~ReconstructionTool();
  
  virtual void execute( Core::ActionContextHandle context );
	/// ACTIVATE:
	/// Activate a tool: this tool is set as the active tool and hence it should
	/// setup the right mouse tools in the viewers.
	virtual void activate();
  
	/// UPDATE_PROGRESS:
	/// Update the progress bar associated with this layer
	void update_progress( double amount, double progress_start = 0.0f, double progress_amount = 1.0f );	
  
	/// UPDATE_PROGRESS_SIGNAL:
	/// When new information on progress is available this signal is triggered. If this signal is 
	/// triggered it should end with a value 1.0 indicating that progress reporting has finished.
	/// Progress is measured between 0.0 and 1.0.
	typedef boost::signals2::signal< void (double) > update_progress_signal_type;
	update_progress_signal_type update_progress_signal_;
  
	/// Whether to show the progress bar
	Core::StateBoolHandle show_progress_bar_state_;
  
	/// Number of iterations
  Core::StateStringHandle input_data_id_;

  Core::StateStringVectorHandle initialGuessSet_state_;

	Core::StateRangedIntHandle iterations_state_;
	Core::StateRangedDoubleHandle measurementScale_state_;
  
	Core::StateStringHandle outputDirectory_state_;	

private:
  ReconstructionToolPrivateHandle private_;
};
  
}

#endif