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

#ifndef APPLICATION_TOOL_ACTIONS_ACTIONNEWMASKLAYER_H
#define APPLICATION_TOOL_ACTIONS_ACTIONNEWMASKLAYER_H


#include <Application/Action/Actions.h>
#include <Application/Interface/Interface.h>
#include <Application/Layer/LayerGroup.h>

namespace Seg3D
{

class ActionNewMaskLayer : public Action
{
	CORE_ACTION( "NewMaskLayer", "<group_name>" );
	
	// -- Constructor/Destructor --
public:
	ActionNewMaskLayer() :
		group_name_("")
	{
		add_parameter("group", group_name_);
	}
	
	virtual ~ActionNewMaskLayer()
	{
	}
	
// -- Functions that describe action --
public:
	virtual bool validate( ActionContextHandle& context );
	virtual bool run( ActionContextHandle& context, ActionResultHandle& result );
	
private:
	ActionParameter< std::string > group_name_;

	// -- Dispatch this action from the interface --
public:
	// CREATE:
	// Create action that moves the layer above
	static ActionHandle Create( LayerGroupHandle group );
	
	// CREATE:
	// Create action that moves the layer above
	static ActionHandle Create( const std::string& group_name );
	
	// DISPATCH
	// Dispatch action that creates a new mask layer 
	static void Dispatch( LayerGroupHandle group );
	
	// DISPATCH
	// Dispatch action that creates a new mask layer 
	static void Dispatch( const std::string& group_name );
	
private:
	// Layer_handle that is requested
	LayerGroupHandle group_handle_;
	
};
	
} // end namespace Seg3D

#endif