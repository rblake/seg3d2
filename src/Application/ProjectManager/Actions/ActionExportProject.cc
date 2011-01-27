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
#include <Application/ProjectManager/Actions/ActionExportProject.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, ExportProject )

namespace Seg3D
{

bool ActionExportProject::validate( Core::ActionContextHandle& context )
{
	if( !ProjectManager::Instance()->current_project_->
		validate_session_name( this->session_name_ ) )
	{
		return false;
	}
	
	boost::filesystem::path path = complete( boost::filesystem::path( 
		this->export_path_.c_str(), boost::filesystem::native ) );

	if( !boost::filesystem::exists( path ) )
	{
		return false;
	}

	return true;
}

bool ActionExportProject::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	bool success = false;

	std::string message = std::string( "Exporting project: '" ) + this->project_name_
		+ std::string( "'" );

	Core::ActionProgressHandle progress = 
		Core::ActionProgressHandle( new Core::ActionProgress( message ) );

	progress->begin_progress_reporting();

	if( ProjectManager::Instance()->export_project( this->export_path_,
		this->project_name_, this->session_name_ ) )
	{
		success = true;
	}

	progress->end_progress_reporting();

	return success;
}

void ActionExportProject::Dispatch( Core::ActionContextHandle context, 
	const std::string& export_path, const std::string& project_name, 
	const std::string& session_name )
{
	ActionExportProject* action = new ActionExportProject;
	
	action->export_path_ = export_path;
	action->project_name_ = project_name;
	action->session_name_ = session_name;
	
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

} // end namespace Seg3D
