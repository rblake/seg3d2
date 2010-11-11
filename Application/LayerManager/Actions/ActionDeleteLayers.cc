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

// Application includes
#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/LayerManager/LayerUndoBuffer.h>
#include <Application/LayerManager/Actions/ActionDeleteLayers.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, DeleteLayers )

namespace Seg3D
{

bool ActionDeleteLayers::validate( Core::ActionContextHandle& context )
{
	std::string error;
	if ( !( LayerManager::CheckGroupExistance( this->group_id_.value(), error ) ) )
	{
		context->report_error( error );
		return false;
	}

	return true; // validated
}

bool ActionDeleteLayers::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	// NOTE: This action will reset the undo buffers as components needed for the undo may have been
	// deleted.
	// NOTE: Most programs assume when a user confirms a delete that one cannot undo there after.
	LayerUndoBuffer::Instance()->reset_undo_buffer();
	
	LayerGroupHandle group = LayerManager::FindLayerGroup( this->group_id_.value() );
	LayerManager::Instance()->delete_layers( group );
	
	return true;
}

Core::ActionHandle ActionDeleteLayers::Create( LayerGroupHandle group )
{
	ActionDeleteLayers* action = new ActionDeleteLayers;
	action->group_id_.value() = group->get_group_id();

	return Core::ActionHandle( action );
}


void ActionDeleteLayers::Dispatch( Core::ActionContextHandle context, LayerGroupHandle group )
{
	Core::ActionDispatcher::PostAction( Create( group ), context );
}

} // end namespace Seg3D
