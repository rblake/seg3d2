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


#include <Application/BackscatterReconstruction/CalibrationTool.h>
#include <Application/BackscatterReconstruction/Actions/ActionCalibrationSegment.h>
#include <Application/BackscatterReconstruction/Actions/ActionCalibrationFitGeometry.h>

#include <Application/Tool/ToolFactory.h>

#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/Layer/LayerManager.h>
#include <Application/ViewerManager/ViewerManager.h>
#include <Application/ViewerManager/Actions/ActionPickPoint.h>

#include <Core/Geometry/Point.h>
#include <Core/Viewer/Mouse.h>
#include <Core/State/Actions/ActionSetAt.h>

#include <Core/Volume/DataVolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

#include <vector>

// test
#include <iostream>
// test

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, CalibrationTool )

namespace Seg3D
{

class CalibrationToolPrivate : public Core::RecursiveLockable
{
public:
  CalibrationTool* tool_;
	void handle_layer_group_insert( LayerHandle layerHandle, bool newGroup );

	ViewerHandle viewer_;
};

void CalibrationToolPrivate::handle_layer_group_insert( LayerHandle layerHandle, bool newGroup )
{
  std::ostringstream oss;
  oss << "layer handle=" << layerHandle->get_layer_name() << ", " << layerHandle->get_layer_id() << ", new group created=" << newGroup;
  CORE_LOG_MESSAGE(oss.str());

  if (layerHandle->get_type() == Core::VolumeType::DATA_E)
  {
    this->tool_->target_layer_state_->set(layerHandle->get_layer_id());
    Core::Point p; // (0, 0, 0) by default
    ActionPickPoint::Dispatch( Core::Interface::GetWidgetActionContext(), p );
  }
  else if (layerHandle->get_type() == Core::VolumeType::MASK_E)
  {
    if (newGroup)
    {
      CORE_LOG_WARNING("Inserting layers from new group");
      return;
    }
    
    this->tool_->maskLayers_state_->add(layerHandle->get_layer_id());
  }
}


CalibrationTool::CalibrationTool( const std::string& toolid ) :
  SingleTargetTool( Core::VolumeType::DATA_E, toolid ),
  private_( new CalibrationToolPrivate )
{
  this->private_->tool_ = this;

  this->add_state( "calibrationSet", this->calibrationSet_state_ );
  // should have 3 layer entries and be found in maskLayers
  this->calibrationSet_state_->add("<none>");
  this->calibrationSet_state_->add("<none>");
  this->calibrationSet_state_->add("<none>");

  this->add_state( "maskLayers", this->maskLayers_state_ );
  
	this->add_connection( LayerManager::Instance()->layer_inserted_signal_.connect( 
    boost::bind( &CalibrationToolPrivate::handle_layer_group_insert, this->private_, _1, _2 ) ) );
}

CalibrationTool::~CalibrationTool()
{
  this->disconnect_all();
}

void CalibrationTool::save( Core::ActionContextHandle context, int index, std::string layerid )
{
  std::cerr << "CalibrationTool::save: " << index << ", " << layerid << std::endl;
  Core::ActionSetAt::Dispatch( context, this->calibrationSet_state_, index, layerid );
}
  
  
void CalibrationTool::execute( Core::ActionContextHandle context )
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

  ActionCalibrationFitGeometry::Dispatch( context,
                                          this->target_layer_state_->get(),
                                          this->calibrationSet_state_->get() );	
}

void CalibrationTool::segment( Core::ActionContextHandle context )
{
  ActionCalibrationSegment::Dispatch( context, this->target_layer_state_->get() );	
}

void CalibrationTool::activate()
{
  Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(),
                             ViewerManager::Instance()->layout_state_,
                             ViewerManager::VIEW_SINGLE_C );
  ViewerHandle viewer = ViewerManager::Instance()->get_active_viewer();
  viewer->view_mode_state_->set(Viewer::AXIAL_C);
}

}