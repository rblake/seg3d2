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

#ifndef APPLICATION_LAYER_MASKLAYER_H
#define APPLICATION_LAYER_MASKLAYER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

// Core includes
#include <Core/Isosurface/Isosurface.h>
#include <Core/Volume/MaskVolume.h>
#include <Core/Utils/AtomicCounter.h>

// Application includes
#include <Application/Layer/Layer.h>

namespace Seg3D
{

// CLASS MaskLayer

// Class definition
class MaskLayer : public Layer
{

	// -- constructor/destructor --
public:

	MaskLayer( const std::string& name, const Core::MaskVolumeHandle& volume );
	MaskLayer( const std::string& state_id );

	// NOTE: This needs to be virtual so it can be deleted from the base class
	virtual ~MaskLayer();

public:
	// TYPE:
	// This function returns the a Core::VolumeType indicating that this is a MaskLayer.
	virtual Core::VolumeType type() const { return Core::VolumeType::MASK_E; }

	// GET_GRID_TRANSFORM:
	// This function returns the grid transform of the mask volume.
	virtual Core::GridTransform get_grid_transform() const;

	// GET_DATA_TYPE:
	// Get the data type of the underlying data
	virtual Core::DataType get_data_type() const;

	// IS_VALID:
	// Check whether the layer has valid data.
	virtual bool is_valid() const;
	
	// GET_VOLUME:
	// Get the underlying volume.
	virtual Core::VolumeHandle get_volume() const;

	// GET_MASK_VOLUME:
	// This function returns the mask volume.
	Core::MaskVolumeHandle get_mask_volume() const;
	
	// SET_MASK_VOLUME:
	// This function set the mask volume to a new mask.
	void set_mask_volume( Core::MaskVolumeHandle volume );

	// GET_ISOSURFACE:
	// Get the iso-surface associated with this layer
	Core::IsosurfaceHandle get_isosurface();

	// COMPUTE_ISOSURFACE
	// Compute the isosurface for this layer using the given quality factor.
	// Quality factor must be one of: 1.0, 0.5, 0.25, 0.125
	void compute_isosurface( double quality_factor );

	void delete_isosurface();

	// -- Signals --
public:
	boost::signals2::signal< void () > isosurface_updated_signal_;
	
	
	// -- state variables --
public:

	// Which color to use for displaying the mask
	Core::StateIntHandle color_state_;

	// State that describes whether to use a fat border for the mask
	Core::StateOptionHandle border_state_;

	// State that describes whether to fill the mask solid
	Core::StateOptionHandle fill_state_;

	// State that describes whether to show the  isosurface state
	Core::StateBoolHandle show_isosurface_state_;
	
	// State that describes whether the iso surface has been generated
	Core::StateBoolHandle iso_generated_state_;

	// -- GUI related states --
public:
	// Whether to display mask visual property control
	Core::StateBoolHandle show_mask_property_state_;

protected:
	// PRE_SAVE_STATES:
	// this function synchronize the generation number for the session saving
	virtual bool pre_save_states( Core::StateIO& state_io );

	// POST_LOAD_STATES:
	// this function takes care of connecting the MaskVolume to the mask layer after it's settings
	// have been loaded
	virtual bool post_load_states( const Core::StateIO& state_io );

	// CLEAN_UP:
	// this function cleans up the mask volume for when you are deleting the mask and reloading 
	virtual void clean_up();

	// -- internal functions --
private:
	void initialize_states();
	void handle_mask_data_changed();
	void handle_isosurface_update_progress( double progress );
	
private:
	// bool to enable a different behavior when the iso surface is being generated after a session
	// load.
	bool loading_;
	
	// Extra private state information
	// NOTE: This used for saving the bit that is used in a mask to a session file. As the state
	// variables are read first, this will allow for reconstructing which data block and which bit
	// need to be loaded.
	Core::StateIntHandle   bit_state_;
	
	// Information about two components not included in the state manager.
	Core::MaskVolumeHandle mask_volume_;
	Core::IsosurfaceHandle isosurface_;
};

} // end namespace Seg3D

#endif
