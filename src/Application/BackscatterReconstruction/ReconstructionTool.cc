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


#include <Application/BackscatterReconstruction/ReconstructionTool.h>
#include <Application/Tool/ToolFactory.h>

#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/Layer/LayerManager.h>
#include <Application/ViewerManager/ViewerManager.h>

#include <Core/Viewer/Mouse.h>

#include <Core/Volume/DataVolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

// test
#include <iostream>
#include <boost/thread.hpp>
// test

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ReconstructionTool )

namespace Seg3D
{

class ReconstructionToolPrivate
{
public:
  void handleOutputDirChanged();

  ReconstructionTool* tool_;
};

void ReconstructionToolPrivate::handleOutputDirChanged() 
{
  std::cerr << "handleOutputDirChanged(): " << tool_->outputDirectory_state_->export_to_string() << std::endl;
}

ReconstructionTool::ReconstructionTool( const std::string& toolid ) :
  SingleTargetTool( Core::VolumeType::MASK_E, toolid ),
  private_( new ReconstructionToolPrivate )
{
	this->private_->tool_ = this;

//  // Create an empty list of label options
//  std::vector< LayerIDNamePair > empty_list( 1, std::make_pair( Tool::NONE_OPTION_C, Tool::NONE_OPTION_C ) );
//  
//  this->add_state( "input_b", this->input_b_state_, Tool::NONE_OPTION_C, empty_list );
//  this->add_extra_layer_input( this->input_b_state_, Core::VolumeType::MASK_E );
//  this->add_state( "input_c", this->input_c_state_, Tool::NONE_OPTION_C, empty_list );
//  this->add_extra_layer_input( this->input_c_state_, Core::VolumeType::MASK_E );
//  this->add_state( "input_d", this->input_d_state_, Tool::NONE_OPTION_C, empty_list );
//  this->add_extra_layer_input( this->input_d_state_, Core::VolumeType::MASK_E );
	this->add_state( "outputDirectory", this->outputDirectory_state_, "" );
  
	// add number of iterations
	this->add_state( "iterations", this->iterations_state_, 20, 1, 100, 1 );
  this->add_connection( this->outputDirectory_state_->state_changed_signal_.connect(
    boost::bind( &ReconstructionToolPrivate::handleOutputDirChanged, this->private_ ) ) );
}

ReconstructionTool::~ReconstructionTool()
{
  this->disconnect_all();
}

void ReconstructionTool::update_progress( double amount, double progress_start, double progress_amount )
{
  this->update_progress_signal_( progress_start + amount * progress_amount );
}
  

void ReconstructionTool::execute( Core::ActionContextHandle context )
{
  // action for saving data goes here...
  for (int i = 0; i < 100; ++i)
  {
    this->update_progress(static_cast<double>(i));
    boost::this_thread::sleep(boost::posix_time::milliseconds(100));
  }
}

void ReconstructionTool::activate()
{
  Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(),
                            ViewerManager::Instance()->layout_state_, ViewerManager::VIEW_SINGLE_C );
  ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
  viewer->view_mode_state_->set(Viewer::AXIAL_C);
}
  
}