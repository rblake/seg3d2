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
#include <Application/BackscatterReconstruction/Actions/ActionReconstructionTool.h>

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
	void handle_layer_group_insert( LayerHandle layerHandle, bool newGroup );

  ReconstructionTool* tool_;
};
  
void ReconstructionToolPrivate::handle_layer_group_insert( LayerHandle layerHandle, bool newGroup )
{
  if (layerHandle->get_type() == Core::VolumeType::DATA_E)
  {
    this->tool_->input_data_id_->set(layerHandle->get_layer_id());
  }
}
  
void ReconstructionToolPrivate::handleOutputDirChanged() 
{
  std::string dir = tool_->outputDirectory_state_->export_to_string();
  if (dir == "[[]]" || dir == "[]") tool_->outputDirectory_state_->set("");
}

ReconstructionTool::ReconstructionTool( const std::string& toolid ) :
  SingleTargetTool( Core::VolumeType::MASK_E, toolid ),
  private_( new ReconstructionToolPrivate )
{
	this->private_->tool_ = this;

	// add number of iterations
  this->add_state( "input_data_id", this->input_data_id_, "<none>" );
	this->add_state( "iterations", this->iterations_state_, 3, 1, 100, 1 );
	this->add_state( "measurementScale", this->measurementScale_state_, 5.0, 1.0, 10.0, 0.1 );
  this->add_state( "initialGuessSet", this->initialGuessSet_state_ );

	this->add_state( "outputDirectory", this->outputDirectory_state_, "" );

  this->add_connection( this->outputDirectory_state_->state_changed_signal_.connect(
    boost::bind( &ReconstructionToolPrivate::handleOutputDirChanged, this->private_ ) ) );
	this->add_connection( LayerManager::Instance()->layer_inserted_signal_.connect( 
    boost::bind( &ReconstructionToolPrivate::handle_layer_group_insert, this->private_, _1, _2 ) ) );
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
  if (this->input_data_id_->get() == "<none>")
  {
    LayerHandle activeLayer = LayerManager::Instance()->get_active_layer();
    if (! activeLayer)
    {
      CORE_LOG_DEBUG("No active layer");
      return;
    }
    LayerGroupHandle groupHandle = activeLayer->get_layer_group();
    if (! groupHandle->has_a_valid_layer() )
    {
      CORE_LOG_DEBUG("Could not find a valid layer in this group");
      return;
    }
    std::vector< LayerHandle > layers;
    groupHandle->get_layers( layers );
    for (size_t i = 0; i < layers.size(); ++i)
    {
      if ( layers[i]->get_type() == Core::VolumeType::DATA_E )
      {
        this->input_data_id_->set( layers[i]->get_layer_id() );
        break;
      }
    }
  }
     
  ActionReconstructionTool::Dispatch( context,
                                     this->input_data_id_->get(),
                                     this->outputDirectory_state_->get(),
                                     this->initialGuessSet_state_->get(),
                                     iterations_state_->get(),
                                     measurementScale_state_->get() );
  // fake something happening
  for (int i = 0; i < 10; ++i)
  {
    this->update_progress(static_cast<double>(i));
    boost::this_thread::sleep(boost::posix_time::milliseconds(100));
  }
  this->update_progress(static_cast<double>(100));
}

void ReconstructionTool::activate()
{
  Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(),
                             ViewerManager::Instance()->layout_state_,
                             ViewerManager::VIEW_SINGLE_C );
  ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
  viewer->view_mode_state_->set(Viewer::AXIAL_C);
}
  
}