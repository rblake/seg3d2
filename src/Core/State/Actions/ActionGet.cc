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

#include <Core/State/StateEngine.h>
#include <Core/State/Actions/ActionGet.h>
#include <Core/Interface/Interface.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Core, Get )

namespace Core
{

bool ActionGet::validate( ActionContextHandle& context )
{
	// Check whether the state exists

	// NOTE: We use lock() to avoid constructor from throwing an exception
	StateBaseHandle state( state_weak_handle_.lock() );

	// If not the state cannot be retrieved report an error
	if ( !state.get() )
	{
		if ( !( StateEngine::Instance()->get_state( stateid_, state ) ) )
		{
			context->report_error( std::string( "Unknown state variable '" ) + stateid_ + "'" );
			return false;
		}
		state_weak_handle_ = state;
	}

	// The action should be able to run
	return true;
}

bool ActionGet::run( ActionContextHandle& context, ActionResultHandle& result )
{
	// Get the state
	StateBaseHandle state( state_weak_handle_.lock() );

	if ( state.get() )
	{
		// Retrieve the current state
		result = ActionResultHandle( new ActionResult );
		state->export_to_variant( *result );
		return true;
	}

	// Signal action is done
	return false;
}


void ActionGet::Dispatch( ActionContextHandle context, StateBaseHandle& state )
{
	// Create new action
	ActionGet* action = new ActionGet;
	action->stateid_ = state->get_stateid();
	action->state_weak_handle_ = state;

	// Post the a action
	ActionDispatcher::PostAction( ActionHandle( action ), context );
}

} // end namespace Core
