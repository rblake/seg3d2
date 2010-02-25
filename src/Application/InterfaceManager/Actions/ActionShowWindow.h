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

#ifndef APPLICATION_INTERFACEMANAGER_ACTIONS_ACTIONSHOWTOOL_H
#define APPLICATION_INTERFACEMANAGER_ACTIONS_ACTIONSHOWTOOL_H

#include <Application/Action/Actions.h>

namespace Seg3D
{

class ActionShowWindow : public Action
{
SCI_ACTION_TYPE("ShowWindow","ShowWindow <windowid>",INTERFACE_E)

	// -- Constructor/Destructor --
public:
	ActionShowWindow()
	{
		add_argument( windowid_ );
	}

	virtual ~ActionShowWindow()
	{
	}

	// -- Functions that describe action --
	virtual bool validate( ActionContextHandle& context );
	virtual bool run( ActionContextHandle& context, ActionResultHandle& result );

	// -- Action parameters --
private:
	ActionParameter< std::string > windowid_;

	// -- Dispatcher for the GUI --
public:

	// CREATE:
	// Create the action
	static ActionHandle Create( const std::string& windowid );

	// DISPATCH:
	// Create the action and dispatch it
	static void Dispatch( const std::string& windowid );

};

} // end namespace Seg3D

#endif
