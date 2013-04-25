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
SCI_REGISTER_TOOL( Seg3D, CalibrationTool )

namespace Seg3D
{

class CalibrationToolPrivate : public Core::RecursiveLockable
{
public:
  CalibrationTool* calibration_tool_;

	ViewerHandle viewer_;

	int center_x_;
	int center_y_;
};



CalibrationTool::CalibrationTool( const std::string& toolid ) :
  SingleTargetTool( Core::VolumeType::MASK_E, toolid ),
  private_( new CalibrationToolPrivate )
{
	// Create an empty list of label options
	std::vector< LayerIDNamePair > empty_list( 1, std::make_pair( Tool::NONE_OPTION_C, Tool::NONE_OPTION_C ) );

	this->add_state( "input_b", this->input_b_state_, Tool::NONE_OPTION_C, empty_list );
	this->add_extra_layer_input( this->input_b_state_, Core::VolumeType::MASK_E );
	this->add_state( "input_c", this->input_c_state_, Tool::NONE_OPTION_C, empty_list );
	this->add_extra_layer_input( this->input_c_state_, Core::VolumeType::MASK_E );
	this->add_state( "input_d", this->input_d_state_, Tool::NONE_OPTION_C, empty_list );
	this->add_extra_layer_input( this->input_d_state_, Core::VolumeType::MASK_E );
}

CalibrationTool::~CalibrationTool()
{
	//this->disconnect_all();
}

bool CalibrationTool::handle_mouse_press( ViewerHandle viewer,  const Core::MouseHistory& mouse_history,
                                         int button, int buttons, int modifiers )
{
	//if ( this->private_->data_layer_ ) return false;
  
	// NOTE: This function call is running on the interface thread
	// Hence we need to lock as most of the paint tool runs on the Application thread.
	{
		CalibrationToolPrivate::lock_type lock( this->private_->get_mutex() );
    
		this->private_->viewer_ = viewer;
   
    // If it is a volume view we cannot do any painting
    if ( this->private_->viewer_->is_volume_view() )
    {
      viewer->set_cursor( Core::CursorShape::ARROW_E );
      return false;
    }
    
    const std::string& target_layer_id = this->target_layer_state_->get();
    Core::VolumeSliceHandle volumeSliceHandle = viewer->get_volume_slice( target_layer_id );
    
    if (! volumeSliceHandle)
    {
      return false;
    }
    
    if ( volumeSliceHandle->volume_type() != Core::VolumeType::MASK_E )
    {
      return false;
    }
    
    Core::MaskVolumeSliceHandle targetSliceHandle = boost::dynamic_pointer_cast< Core::MaskVolumeSlice >( volumeSliceHandle );
    if ( ! targetSliceHandle->is_valid() )
    {
      return false;
    }
    
		this->private_->viewer_ = viewer;
    
		double xpos, ypos;
		viewer->window_to_world( mouse_history.current_.x_, mouse_history.current_.y_, xpos, ypos );
    
		int i, j;
		targetSliceHandle->world_to_index( xpos, ypos, i, j );    
	}

	// Pass handling on to normal handler 
  return false;
}
  
bool CalibrationTool::handle_mouse_move( ViewerHandle viewer, const Core::MouseHistory& mouse_history,
                       int button, int buttons, int modifiers )
{
	//if ( this->private_->data_layer_ ) return false;
  
	// NOTE: This function call is running on the interface thread
	// Hence we need to lock as most of the tool runs on the Application thread.
	{
		CalibrationToolPrivate::lock_type lock( this->private_->get_mutex() );
  
    // If it is a volume view we cannot use the tool
    if ( viewer->is_volume_view() )
    {
      viewer->set_cursor( Core::CursorShape::ARROW_E );
      return false;
    }

//    if( !( modifiers == Core::KeyModifier::NO_MODIFIER_E || 
//          modifiers == Core::KeyModifier::CONTROL_MODIFIER_E ||
//          modifiers == Core::KeyModifier::GROUPSWITCH_MODIFIER_E ) )
//    {
//      return false;
//    }
//
//    if( button != Core::MouseButton::LEFT_BUTTON_E )
//    {
//      return false;
//    }

    const std::string& target_layer_id = this->target_layer_state_->get();
    std::cerr << "target_layer_id=" << target_layer_id << std::endl;
    Core::VolumeSliceHandle volumeSliceHandle = viewer->get_volume_slice( target_layer_id );
    
    if (! volumeSliceHandle)
    {
      return false;
    }
    
    if ( volumeSliceHandle->volume_type() != Core::VolumeType::MASK_E )
    {
      return false;
    }
    
    Core::MaskVolumeSliceHandle targetSliceHandle = boost::dynamic_pointer_cast< Core::MaskVolumeSlice >( volumeSliceHandle );
    if ( ! targetSliceHandle->is_valid() )
    {
      return false;
    }

		this->private_->viewer_ = viewer;

		double xpos, ypos;
		viewer->window_to_world( mouse_history.current_.x_, mouse_history.current_.y_, xpos, ypos );
		int i, j;
		targetSliceHandle->world_to_index( xpos, ypos, i, j );

    bool value = targetSliceHandle->get_mask_at( static_cast<size_t>( i ), static_cast<size_t>( j ) );
    if (value)
    {
      viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );
    }
    else
    {
      viewer->set_cursor( Core::CursorShape::CROSS_E );
    }
	}

	// Pass handling on to normal handler 
  return false;
}
  
void CalibrationTool::execute( Core::ActionContextHandle context )
{
  // action for saving data goes here...
  for (int i = 0; i < 1000; ++i)
  {
    // do nothing
  }
}


}