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

#ifndef APPLICATION_TOOL_ACTIONS_ACTIONCONFIDENCECONNECTEDFILTER_H
#define APPLICATION_TOOL_ACTIONS_ACTIONCONFIDENCECONNECTEDFILTER_H

#include <Core/Action/Actions.h>
#include <Core/Interface/Interface.h>
#include <Application/Layer/Layer.h>

namespace Seg3D
{
	
class ActionConfidenceConnectedFilter : public Core::Action
{

CORE_ACTION_XML(
	CORE_ACTION_TYPE( "ConfidenceConnectedFilter", "Find regions that are connected within a value range." )
	CORE_ACTION_ARGUMENT( "layerid", "The layerid on which this filter needs to be run." )
	CORE_ACTION_KEY( "iterations", "10", "Number of iterations to recomputed the range based on the last selected region." )
	CORE_ACTION_KEY( "multiplier", "1", "Multiplier" )
)
	
	// -- Constructor/Destructor --
public:
	ActionConfidenceConnectedFilter()
	{
		add_argument( this->layer_id_ );
		
		add_key( this->iterations_ );
		add_key( this->multiplier_ );
		
		add_cachedhandle( this->layer_ );
	}
	
	virtual ~ActionConfidenceConnectedFilter()
	{
	}
	
	// -- Functions that describe action --
public:
	virtual bool validate( Core::ActionContextHandle& context );
	virtual bool run( Core::ActionContextHandle& context, Core::ActionResultHandle& result );
	
	// -- Action parameters --
private:

	Core::ActionParameter< std::string > layer_id_;
	Core::ActionParameter< int > iterations_;
	Core::ActionParameter< int > multiplier_;
	
	Core::ActionCachedHandle< LayerHandle > layer_;
	
	// -- Dispatch this action from the interface --
public:

	// CREATE:	
	// Create the action, but do not dispatch it
	static Core::ActionHandle Create( std::string layer_id, int iterations, int multiplier );
			
	// DISPATCH:
	// Create and dispatch action that inserts the new layer 
	static void Dispatch( Core::ActionContextHandle context, std::string layer_id, 
		int iterations, int multiplier );
	
};
	
} // end namespace Seg3D

#endif