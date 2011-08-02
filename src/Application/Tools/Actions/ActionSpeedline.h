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

#ifndef APPLICATION_TOOLS_ACTIONS_ACTIOINSPEEDLINE_H
#define APPLICATION_TOOLS_ACTIONS_ACTIOINSPEEDLINE_H

// Core includes
#include <Core/Volume/VolumeSlice.h>
#include <Core/Utils/AtomicCounter.h>

// Application includes
#include <Application/Layer/LayerAction.h>

namespace Seg3D
{

class ActionSpeedlinePrivate;
typedef boost::shared_ptr< ActionSpeedlinePrivate > ActionSpeedlinePrivateHandle;

class ActionSpeedline : public LayerAction
{

CORE_ACTION
( 
	CORE_ACTION_TYPE( "Speedline", "Fill or erase a slice of a mask layer within "
										"the region enclosed by the Speedline.")
	CORE_ACTION_ARGUMENT( "target", "The ID of the target data layer." )
	CORE_ACTION_ARGUMENT( "slice_type", "The slicing direction to be painted on." )
	CORE_ACTION_ARGUMENT( "slice_number", "The slice number to be painted on." )
	CORE_ACTION_ARGUMENT( "vertices", "The 2D coordinates of Speedline vertices." )
	CORE_ACTION_ARGUMENT( "current_vertex_index", "The vertex needes to compute paths." )
	CORE_ACTION_OPTIONAL_ARGUMENT( "iterations", "1000", "Number of iterations to perform." )
	CORE_ACTION_OPTIONAL_ARGUMENT( "termination", "1.0", "Unit of Termination." )
	CORE_ACTION_OPTIONAL_ARGUMENT( "update_all_path", "true", "Update all paths" )
	CORE_ACTION_OPTIONAL_ARGUMENT( "itk_path_state_id", "", "The statid of the state variable into which ITK continuous index values will be written." )
	CORE_ACTION_OPTIONAL_ARGUMENT( "world_path_state_id", "", "The statid of the state variable into which world coordinate path values will be written." )
	CORE_ACTION_OPTIONAL_ARGUMENT( "path_vertices_state_id", "", "The stateid of the state variable into which vertices values will be written." )
	CORE_ACTION_CHANGES_PROJECT_DATA()
)

public:
	ActionSpeedline();

	// VALIDATE:
	// Each action needs to be validated just before it is posted. This way we
	// enforce that every action that hits the main post_action signal will be
	// a valid action to execute.
	virtual bool validate( Core::ActionContextHandle& context );

	// RUN:
	// Each action needs to have this piece implemented. It spells out how the
	// action is run. It returns whether the action was successful or not.
	virtual bool run( Core::ActionContextHandle& context, Core::ActionResultHandle& result );

private:
	ActionSpeedlinePrivateHandle private_;

public:

	static void Dispatch( Core::ActionContextHandle context, const std::string& layer_id,
		Core::VolumeSliceType slice_type,
		size_t slice_number,
		const std::vector< Core::Point > vertices, 
		int current_vertex_index,
		int iterations, 
		double termination,
		bool update_all_paths,
		const std::string& itk_path_state_id,
		const std::string& world_path_state_id,
		const std::string& path_vertices_state_id,
		long  action_id,
		Core::AtomicCounterHandle action_handle
		);
};

} // end namespace Seg3D

#endif
