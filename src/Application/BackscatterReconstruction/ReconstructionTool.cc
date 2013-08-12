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

#include <Core/State/Actions/ActionSet.h>

//#include <Core/Viewer/Mouse.h>

#include <Core/Volume/DataVolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

// test
#include <iostream>
// test

#include <boost/filesystem.hpp>

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ReconstructionTool )

namespace Seg3D
{

class ReconstructionToolPrivate
{
public:
  ReconstructionToolPrivate() {}
  ~ReconstructionToolPrivate() {}

  void handleOutputDirChanged();
	void handle_layer_group_insert( LayerHandle layerHandle, bool newGroup );
	void handle_abort();
  
	void update_progress( double amount, double progress_start = 0.0f, double progress_amount = 1.0f );
  void reset_progress();

  ReconstructionTool* tool_;
};
  
void ReconstructionToolPrivate::handle_layer_group_insert( LayerHandle layerHandle, bool newGroup )
{
  if (layerHandle->get_type() == Core::VolumeType::DATA_E)
  {
    this->tool_->target_layer_state_->set(layerHandle->get_layer_id());
  }
  else if (layerHandle->get_type() == Core::VolumeType::MASK_E)
  {
    if ( this->tool_->input_b_state_->get() == "<none>" )
    {
      this->tool_->input_b_state_->set(layerHandle->get_layer_id());
    }
    else if ( this->tool_->input_c_state_->get() == "<none>" )
    {
      this->tool_->input_c_state_->set(layerHandle->get_layer_id());
    }
  }
}

void ReconstructionToolPrivate::handle_abort()
{
  CORE_LOG_MESSAGE("Abort called from reconstruction tool");
  ActionReconstructionTool::Abort();
}

void ReconstructionToolPrivate::handleOutputDirChanged() 
{
  std::string dir = tool_->outputDirectory_state_->export_to_string();
  if (dir == "[[]]" || dir == "[]") tool_->outputDirectory_state_->set("");
}

ReconstructionTool::ReconstructionTool( const std::string& toolid ) :
  SingleTargetTool( Core::VolumeType::DATA_E/*|Core::VolumeType::MASK_E*/, toolid ),
  private_( new ReconstructionToolPrivate )
{
	this->private_->tool_ = this;

	// Create an empty list of label options
	std::vector< LayerIDNamePair > empty_list( 1, 
    std::make_pair( Tool::NONE_OPTION_C, Tool::NONE_OPTION_C ) );
  
	// add number of iterations
	this->add_state( "iterations", this->iterations_state_, 1, 0, 3, 1 );
	this->add_state( "xyVoxelSizeScale", this->xyVoxelSizeScale_state_, 0.5, 0.1, 10.0, 0.1 );
	this->add_state( "zVoxelSizeScale", this->zVoxelSizeScale_state_, 0.5, 0.1, 10.0, 0.1 );
	this->add_state( "outputDirectory", this->outputDirectory_state_, "" );

	this->add_state( "input_b", this->input_b_state_, Tool::NONE_OPTION_C, empty_list );
	this->add_extra_layer_input( this->input_b_state_, Core::VolumeType::MASK_E, false, false );
	this->add_state( "input_c", this->input_c_state_, Tool::NONE_OPTION_C, empty_list );
	this->add_extra_layer_input( this->input_c_state_, Core::VolumeType::MASK_E, false, false );
  
  this->add_connection( this->outputDirectory_state_->state_changed_signal_.connect(
    boost::bind( &ReconstructionToolPrivate::handleOutputDirChanged, this->private_ ) ) );

	this->add_connection( LayerManager::Instance()->layer_inserted_signal_.connect( 
    boost::bind( &ReconstructionToolPrivate::handle_layer_group_insert, this->private_, _1, _2 ) ) );

	this->add_connection( this->abort_signal_.connect(
    boost::bind( &ReconstructionToolPrivate::handle_abort, this->private_ ) ) );
}

void ReconstructionToolPrivate::update_progress( double amount, double progress_start, double progress_amount )
{
  std::ostringstream oss;
  oss << "Update reconstruction progress: " << amount;
  CORE_LOG_DEBUG(oss.str());
  this->tool_->update_progress_signal_( progress_start + amount * progress_amount );
}
  
void ReconstructionToolPrivate::reset_progress()
{
  CORE_LOG_DEBUG("Reset reconstruction progress");
  this->tool_->reset_progress_signal_();
}

ReconstructionTool::~ReconstructionTool()
{
  this->disconnect_all();
}
  

void ReconstructionTool::execute( Core::ActionContextHandle context )
{
  if ( this->target_layer_state_->get() == Tool::NONE_OPTION_C )
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
        this->target_layer_state_->set( layers[i]->get_layer_id() );
        break;
      }
    }
  }

  std::vector< std::string > initialGuessSet;
  if ( this->input_b_state_->get() != Tool::NONE_OPTION_C )
  {
    initialGuessSet.push_back(this->input_b_state_->get());
  }

  if ( this->input_c_state_->get() != Tool::NONE_OPTION_C )
  {
    initialGuessSet.push_back(this->input_c_state_->get());
  }
  
  if (initialGuessSet.size() > 0 && initialGuessSet.size() < 2)
  {
    CORE_LOG_WARNING("Need 2 segmentations for the initial reconstruction. Running reconstruction without initial reconstruction data.");
    initialGuessSet.clear();
  }

  if ( this->outputDirectory_state_->get().size() == 0 || ! boost::filesystem::is_directory( this->outputDirectory_state_->get() ) )
  {
    boost::filesystem::path algorithm_work_dir;
    boost::filesystem::path algorithm_config_file;    
    boost::filesystem::path algorithm_source_illum_file;
    boost::filesystem::path algorithm_geometry_file;
    Core::Application::Instance()->get_algorithm_config(algorithm_work_dir,
                                                        algorithm_config_file,
                                                        algorithm_source_illum_file,
                                                        algorithm_geometry_file);

    Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(),
                              this->outputDirectory_state_,
                              algorithm_work_dir.string() );
  }

  this->private_->reset_progress();

  ActionReconstructionTool::Dispatch( context,
    this->target_layer_state_->get(),
    initialGuessSet,
    this->outputDirectory_state_->get(),
    iterations_state_->get(),
    xyVoxelSizeScale_state_->get(),
    zVoxelSizeScale_state_->get(),
    boost::bind( &ReconstructionToolPrivate::update_progress, this->private_, _1, _2, _3 ) );
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