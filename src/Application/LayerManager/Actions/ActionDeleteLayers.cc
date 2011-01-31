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
#include <Application/UndoBuffer/UndoBuffer.h>
#include <Application/LayerManager/Actions/ActionDeleteLayers.h>
#include <Application/LayerManager/LayerUndoBufferItem.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, DeleteLayers )

namespace Seg3D
{

bool ActionDeleteLayers::validate( Core::ActionContextHandle& context )
{
	if( this->layers_ != "" )
	{
		std::vector< std::string > layer_vector;
		layer_vector = Core::SplitString( this->layers_, "|" );
	
		for ( size_t j = 0; j < layer_vector.size(); j++ )
		{	
			if ( !( LayerManager::CheckLayerExistance( layer_vector[ j ], 
				context ) ) ) return false;
		}
		
		if( layer_vector.size() ) return true;
	}
	
	context->report_error( "No valid layer ids can be found in the layerid input." );
	return false; 
}

bool ActionDeleteLayers::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	std::vector< std::string > layer_vector;
	layer_vector = Core::SplitString( this->layers_, "|" );
		
	// Create an undo item for this action
	LayerUndoBufferItemHandle item( new LayerUndoBufferItem( "Delete layer(s)" ) );
	// Tell which action has to be re-executed to obtain the result
	item->set_redo_action( this->shared_from_this() );
	// Tell which layers are to be deleted so they can be added back 
	for ( size_t i = 0; i < layer_vector.size(); ++i )
	{
		LayerHandle layer = LayerManager::Instance()->get_layer_by_id( layer_vector[ i ] );
		layer->abort_signal_();
		item->add_layer_to_add( layer );
	}
	// Add the complete record to the undo buffer
	UndoBuffer::Instance()->insert_undo_item( context, item );

	LayerManager::Instance()->delete_layers( layer_vector );
	
	return true;
}

void ActionDeleteLayers::Dispatch( Core::ActionContextHandle context, std::vector< std::string > layers )
{
	ActionDeleteLayers* action = new ActionDeleteLayers;
	action->layers_ = Core::ExportToString( layers );
	
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

} // end namespace Seg3D
