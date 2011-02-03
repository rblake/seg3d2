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
#include <Application/ProjectManager/Actions/ActionSaveProjectAs.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, SaveProjectAs )

namespace Seg3D
{

bool ActionSaveProjectAs::validate( Core::ActionContextHandle& context )
{
	boost::filesystem::path path = complete( boost::filesystem::path( 
		this->export_path_.c_str(), boost::filesystem::native ) );
		
	boost::filesystem::path current_project_path = boost::filesystem::path( 
		ProjectManager::Instance()->current_project_path_state_->get() ) /
		ProjectManager::Instance()->current_project_->project_name_state_->get() /
		( ProjectManager::Instance()->current_project_->project_name_state_->get() + ".s3d" );
		
	if( !boost::filesystem::exists( current_project_path ) && ProjectManager::Instance()->project_saved_state_->get() )
	{
		CORE_LOG_ERROR( "Your project could not be found at the place it was last saved!" );
		ProjectManager::Instance()->project_saved_state_->set( false );
	}

	if( !boost::filesystem::exists( path ) )
	{
		CORE_LOG_ERROR( "Project saving FAILED! The location specified does not exist." );
		return false;
	}

	return true;
}

bool ActionSaveProjectAs::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	bool success = false;

	std::string message = std::string( "Saving project as: '" ) + this->project_name_
		+ std::string( "'" );

	Core::ActionProgressHandle progress = 
		Core::ActionProgressHandle( new Core::ActionProgress( message ) );

	progress->begin_progress_reporting();
	
	boost::filesystem::path path = complete( boost::filesystem::path( 
		this->export_path_.c_str(), boost::filesystem::native ) );

	if( ProjectManager::Instance()->project_save_as( path,
		this->project_name_ ) )
	{
		success = true;
	}

	progress->end_progress_reporting();
	
	if( success )
	{
		ProjectManager::Instance()->get_current_project()->reset_project_changed();
	}
	
	return success;
}

void ActionSaveProjectAs::Dispatch( Core::ActionContextHandle context, 
	const std::string& export_path, const std::string& project_name )
{
	ActionSaveProjectAs* action = new ActionSaveProjectAs;
	
	action->export_path_ = export_path;
	action->project_name_ = project_name;
	
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

} // end namespace Seg3D
