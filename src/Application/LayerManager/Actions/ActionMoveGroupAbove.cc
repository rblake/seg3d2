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


#include <Application/LayerManager/LayerManager.h>
#include <Application/LayerManager/Actions/ActionMoveGroupAbove.h>

namespace Seg3D
{
	
	// REGISTER ACTION:
	// Define a function that registers the action. The action also needs to be
	// registered in the CMake file.
	CORE_REGISTER_ACTION( MoveGroupAbove );
	
	bool ActionMoveGroupAbove::validate( ActionContextHandle& context )
	{

		if ( !( StateEngine::Instance()->is_stateid( this->group_to_move_id_.value() ) ) )
		{
			context->report_error( std::string( "GroupID '" ) + this->group_to_move_id_.value() 
				+ "' is invalid" );
			return false;
		}
		
		if ( !( StateEngine::Instance()->is_stateid( this->group_below_id_.value() ) ) )
		{
			context->report_error( std::string( "GroupID '" ) + this->group_below_id_.value()
				+ "' is invalid" );
			return false;
		}

		return true;
	}
	
	bool ActionMoveGroupAbove::run( ActionContextHandle& context, ActionResultHandle& result )
	{
		if ( ( StateEngine::Instance()->is_stateid( this->group_below_id_.value() ) ) && 
			( StateEngine::Instance()->is_stateid( this->group_to_move_id_.value() ) ) )
		{
			return LayerManager::Instance()->move_group_above( this->group_to_move_id_.value(),
				this->group_below_id_.value() );
		}

		return false;
	}
	
	ActionHandle ActionMoveGroupAbove::Create( const std::string& group_to_move_id, 
		const std::string& group_below_id )
	{
		// Create new action
		ActionMoveGroupAbove* action = new ActionMoveGroupAbove;
		
		// We need to fill in these to ensure the action can be replayed
		action->group_below_id_.value() = group_below_id;
		action->group_to_move_id_.value() = group_to_move_id;
		
		// Post the new action
		return ActionHandle( action );
	}
	

	void ActionMoveGroupAbove::Dispatch( const std::string& group_to_move_id, 
		const std::string& group_below_id )
	{
		Interface::PostAction( Create( group_to_move_id, group_below_id ) );
	}
	
} // end namespace Seg3D