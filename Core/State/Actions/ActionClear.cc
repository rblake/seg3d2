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

#include <Core/Action/ActionFactory.h>
#include <Core/State/Actions/ActionClear.h>
#include <Core/State/StateEngine.h>

CORE_REGISTER_ACTION( Core, Clear )

namespace Core
{

ActionClear::ActionClear()
{
	this->add_argument( this->stateid_ );
}

ActionClear::~ActionClear()
{
}

bool ActionClear::validate( ActionContextHandle& context )
{
	StateBaseHandle state( this->state_weak_handle_.lock() );
	if ( !state )
	{
		if ( !StateEngine::Instance()->get_state( this->stateid_.value(), state ) )
		{
			context->report_error( std::string( "Unknown state variable '" ) + stateid_.value() + "'" );
			return false;
		}
		StateVectorBaseHandle vector_state = 
			boost::dynamic_pointer_cast< StateVectorBase >( state );
		if ( !vector_state )
		{
			context->report_error( std::string( "State variable '") + this->stateid_.value() +
				"' doesn't support ActionClear" );
			return false;
		}
		this->state_weak_handle_ = vector_state;
	}
	
	return true;
}

bool ActionClear::run( ActionContextHandle& context, ActionResultHandle& result )
{
	StateVectorBaseHandle vector_state( this->state_weak_handle_.lock() );
	if ( vector_state )
	{
		vector_state->clear( context->source() );
		return true;
	}

	return false;
}

ActionHandle ActionClear::Create( const StateVectorBaseHandle& state )
{
	ActionClear* action = new ActionClear;
	action->stateid_.set_value( state->stateid() );
	action->state_weak_handle_ = state;

	return ActionHandle( action );
}

void ActionClear::Dispatch( ActionContextHandle context, const StateVectorBaseHandle& state )
{
	ActionDispatcher::PostAction( Create( state ), context );
}

} // end namespace Core
