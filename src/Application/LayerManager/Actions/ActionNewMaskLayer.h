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

#include <Core/Action/Actions.h>
#include <Core/Interface/Interface.h>
#include <Application/Layer/LayerGroup.h>

namespace Seg3D
{

class ActionNewMaskLayer : public Core::Action
{
	CORE_ACTION( "NewMaskLayer", "NewMaskLayer <group_name>" );
	
	// -- Constructor/Destructor --
public:
	ActionNewMaskLayer() :
		group_name_("")
	{
		add_parameter("group", group_name_);
		add_cachedhandle( group_handle_ );
	}
	
	virtual ~ActionNewMaskLayer()
	{
	}
	
// -- Functions that describe action --
public:
	virtual bool validate( Core::ActionContextHandle& context );
	virtual bool run( Core::ActionContextHandle& context, 
		Core::ActionResultHandle& result );
	
private:
	// The name of the group where the mask needs to be added
	Core::ActionParameter< std::string > group_name_;
	
	//  A short cut to where the layer group
	Core::ActionCachedHandle<LayerGroupHandle> group_handle_;

	// -- Dispatch this action from the interface --
public:
	// CREATE:
	// Create action that moves the layer above
	static Core::ActionHandle Create( LayerGroupHandle group );
	
	// DISPATCH
	// Dispatch action that creates a new mask layer 
	static void Dispatch( LayerGroupHandle group );
	
};
	
} // end namespace Seg3D

#endif
