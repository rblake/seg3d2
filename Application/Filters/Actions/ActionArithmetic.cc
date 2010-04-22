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
#include <Application/Filters/Actions/ActionArithmetic.h>

namespace Seg3D
{
	
// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Arithmetic );

bool ActionArithmetic::validate( ActionContextHandle& context )
{
	if( !( StateEngine::Instance()->is_statealias( this->layer_a_alias_ ) ) )
	{
		context->report_error( std::string( "LayerID '" ) + this->layer_a_alias_ + "' is invalid" );
		return false;
	}
	
	if( !( StateEngine::Instance()->is_statealias( this->layer_b_alias_ ) ) )
	{
		context->report_error( std::string( "LayerID '" ) + this->layer_b_alias_ + "' is invalid" );
		return false;
	}
	
	if( !( StateEngine::Instance()->is_statealias( this->layer_c_alias_ ) ) )
	{
		context->report_error( std::string( "LayerID '" ) + this->layer_c_alias_ + "' is invalid" );
		return false;
	}
	
	return true;
}

bool ActionArithmetic::run( ActionContextHandle& context, ActionResultHandle& result )
{
	if ( StateEngine::Instance()->is_statealias( this->layer_a_alias_ ) )
	{
		// TODO: run filter
		context->report_message( "The Arithmetic Filter has been triggered "
			"successfully on: "  + this->layer_a_alias_ + ", " + this->layer_b_alias_
			+ ", and " + this->layer_c_alias_ );
		
		return true;
	}
		
	return false;
}


void ActionArithmetic::Dispatch( std::string layer_a_alias, std::string layer_b_alias, 
								std::string layer_c_alias, std::string expression, bool replace )
{
	ActionArithmetic* action = new ActionArithmetic;
	action->layer_a_alias_ = layer_a_alias;
	action->layer_b_alias_ = layer_b_alias;
	action->layer_c_alias_ = layer_c_alias;
	action->expression_ = expression;
	action->replace_ = replace;
	
	Interface::PostAction( ActionHandle( action ) );
}
	
} // end namespace Seg3D
