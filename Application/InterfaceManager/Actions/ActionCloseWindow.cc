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

#include <Application/Interface/Interface.h>
#include <Application/InterfaceManager/InterfaceManager.h>
#include <Application/InterfaceManager/Actions/ActionCloseWindow.h>

namespace Seg3D
{

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( CloseWindow );

bool ActionCloseWindow::validate( ActionContextHandle& context )
{
	if ( !( InterfaceManager::Instance()->is_windowid( windowid_.value() ) ) )
	{
		context->report_error( std::string( "WindowID '" ) + windowid_.value() + "' is invalid" );
		return ( false );
	}

	return ( true ); // validated
}

bool ActionCloseWindow::run( ActionContextHandle& context, ActionResultHandle& result )
{
	InterfaceManager::Instance()->close_window_signal_( windowid_.value() );
	return ( true ); // success
}

ActionHandle ActionCloseWindow::Create( const std::string& windowid )
{
	// Create new action
	ActionCloseWindow* action = new ActionCloseWindow;

	// Set action parameters
	action->windowid_.value() = windowid;

	// Create action handle
	return ( ActionHandle( action ) );
}

void ActionCloseWindow::Dispatch( const std::string& windowid )
{
	// Post the new action
	Interface::PostAction( Create( windowid ) );
}

} // end namespace Seg3D
