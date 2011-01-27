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

// Boost includes
#include <boost/filesystem.hpp>

// Application includes
#include <Application/ProjectManager/ProjectManager.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/UndoBuffer/UndoBuffer.h>
#include <Application/ToolManager/ToolManager.h>
#include <Application/ProjectManager/Actions/ActionQuickOpen.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, QuickOpen )

namespace Seg3D
{

bool ActionQuickOpen::validate( Core::ActionContextHandle& context )
{
	return true;
}

bool ActionQuickOpen::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	ProjectManager::Instance()->new_project( "", "", false );
	if ( ProjectManager::Instance()->get_current_project() )
	{
		ProjectManager::Instance()->get_current_project()->reset_project_changed();
	}
	
	// Clear undo buffer
	UndoBuffer::Instance()->reset_undo_buffer();
	
	return true;
}

void ActionQuickOpen::Dispatch( Core::ActionContextHandle context )
{
	ActionQuickOpen* action = new ActionQuickOpen;
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

} // end namespace Seg3D
