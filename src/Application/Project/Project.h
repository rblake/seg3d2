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

#ifndef APPLICATION_PROJECT_PROJECT_H
#define APPLICATION_PROJECT_PROJECT_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

// STL includes
#include <string>
#include <vector>

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>

// Sqlite includes
#include <Externals/sqlite/sqlite3.h>

// Application includes
#include <Application/Provenance/Provenance.h>
#include <Application/Provenance/ProvenanceStep.h>
#include <Application/Session/Session.h>
#include <Application/Project/DataManager.h>
#include <Application/DatabaseManager/DatabaseManager.h>

// Core includes
#include <Core/Action/Action.h>
#include <Core/Application/Application.h>
#include <Core/Interface/Interface.h>
#include <Core/Volume/Volume.h>
#include <Core/State/State.h>

namespace Seg3D
{

// CLASS Project
// This is the main class for collecting state information on a Project
class Project;

class SessionInfo
{
public:
	SessionInfo( std::string session_name = "", std::string username = "", std::string timestamp = "" ) :
		session_name_( session_name ),
		username_( username ),
		timestamp_( timestamp )
	{
	}

	virtual ~SessionInfo()
	{
	}

public:
	std::string session_name_;
	std::string username_;
	std::string timestamp_;
};	
	
typedef boost::shared_ptr< Project > ProjectHandle;

// Class definition
class Project : public Core::StateHandler, private DatabaseManager // don't need Core::RecursiveLockable because DatabaseManager is one
{

	// -- constructor/destructor --
public:
	Project( const std::string& project_name );
	virtual ~Project();
	
public:
	Core::StateStringHandle project_name_state_;
	Core::StateBoolHandle save_custom_colors_state_;
	Core::StateStringVectorHandle sessions_state_;
	Core::StateStringVectorHandle project_notes_state_;
	Core::StateLongLongHandle project_file_size_state_;
	Core::StateIntHandle session_count_state_;
	Core::StateStringHandle current_session_name_state_;
	std::vector< Core::StateColorHandle > color_states_;
	
	// Generation counter state, this one is filled out when the project is saved
	// So the generation state can be restored 
	Core::StateLongLongHandle generation_count_state_;
	
	
public:
	typedef boost::signals2::signal< void() > sessions_changed_signal_type;
	typedef boost::signals2::signal<  void( std::vector< std::pair< ProvenanceID, std::string > > ) > provenance_records_signal_type;

	sessions_changed_signal_type sessions_changed_signal_;
	provenance_records_signal_type provenance_record_signal_;

public:
	// INITIALIZE_FROM_FILE:
	// this file initializes the state values for project from the file at the path specified
	bool initialize_from_file( const std::string& project_name );

	// LOAD_SESSION:
	// this function will be called to load a specific session
	bool load_session( const std::string& session_name );
	
	// SAVE_SESSION:
	// this function will be called from the project manager to save a session
	bool save_session(  const std::string& timestamp, const std::string& session_name );
	
	// DELETE_SESSION:
	// this function will be called by the project manager to delete a session
	bool delete_session( const std::string& session_name );

	// PROJECT_EXPORT:
	// this function will export the current project and the passed vector of session names to file
	bool project_export( boost::filesystem::path path, const std::string& project_name, 
		const std::string& session_name );
	
	// SAVE_AS:
	// this function will save the current project as a new project
	bool save_as( boost::filesystem::path path, const std::string& project_name );

	// VALIDATE_SESSION_NAME:
	// function for validating that a session name exists
	bool validate_session_name( std::string& session_name );

	// INVALIDATE_CURRENT_SESSION: // NOT CURRENTLY USED //
	// this is a public function that enables the ProjectManager to call invalidate on the current
	// session
	void invalidate_current_session(){ this->current_session_->invalidate(); }

	// CLEAR_DATAMANAGER_LIST:
	// function for clearing out the datamanager list
	//void clear_datamanager_list(){ this->data_manager_->clear_data_file_list(); }

	// SET_PROJECT_PATH:
	// function that lets the project manager set the project path for the project
	void set_project_path( const boost::filesystem::path& project_path );

	// SET_SIGNAL_BLOCK:
	// this function is a public function that enables the project manager to disable the signals 
	// that the project emits when it's state variables are changed
	void set_signal_block( bool on_off );
	
	// CHECK_PROJECT_CHANGED:
	// Check whether the project was changed
	bool check_project_changed();

protected:
	// PRE_SAVE_STATES:
	// this function synchronizes the colors if they are set to be saved with the project
	virtual bool pre_save_states( Core::StateIO& state_io );
	
	// POST_LOAD_STATES:
	// this function sets Seg3d's mask colors if they are set to be saved with the project
	virtual bool post_load_states( const Core::StateIO& state_io );
	
private:
	// CLEANUP_SESSION_DATABASE:
	// this function cleans up sessions in the session list that have been deleted by the user
	void cleanup_session_database();
	
	// IMPORT_OLD_SESSION_INFO_INTO_DATABASE:
	// this function is for backwards compatibility with older versions of seg3d that store the
	// list of sessions in the project xml file rather than the database
	void import_old_session_info_into_database();
	
	// GET_MOST_RECENT_SESSION:
	// this function true or false based on whether it was able to find the most recent
	// session name
	bool get_most_recent_session_name( std::string& session_name );
	
	// GET_DATA_FILE_SIZE:
	// this function calculates the file size of the data files that are part of the project
	long long get_data_file_size();
	
	void delete_unused_data_files();

	// -- provenance support --
public:	
	// ADD_TO_PROVENANCE_DATABASE:
	// adds the provenance step to the database
	bool add_to_provenance_database( ProvenanceStepHandle& step );

	// GET_PROVENANCE_RECORD:
	// returns a vector that is the provenance record for a particular ProvenanceID
	void request_signal_provenance_record( ProvenanceID prov_id );
	
public:
	// SET_PROJECT_CHANGED:
	// Set that the session has been modified
	void set_project_changed( Core::ActionHandle action, Core::ActionResultHandle result );

	// RESET_PROJECT_CHANGED:
	// Reset the flag that remembers that a session has changed
	void reset_project_changed();
	
	// CREATE_DATABASE_SCHEMA:
	// this is an inherited function that
	virtual bool create_database_schema(); 
	
	// INSERT_SESSION_INTO_DATABASE:
	// this inserts a session into the database
	bool insert_session_into_database( const std::string& timestamp, const std::string& session_name, std::string& error_message, const std::string& user_name = "" );
	
	// DELETE_SESSION_FROM_DATABASE:
	// this deletes a session from the database
	bool delete_session_from_database( const std::string& session_name );
	
	// GET_ALL_SESSIONS:
	// this returns all the sessions
	bool get_all_sessions( std::vector< SessionInfo >& sessions );
	
	// GET_SESSION:
	// this function returns true or false based on whether it can find the session name
	// that was passed it
	bool get_session( SessionInfo& session, const std::string& session_name );
	
	// CLOSE_PROVENANCE_DATABASE:
	// this fucntion closes the provenance database
	void close_provenance_database();
	
	// CHECKPOINT_PROVENANCE_DATABASE:
	// this function dumps the database to disk
	void checkpoint_provenance_database();

private:
	// PROVENANCE_RECURSIVE_HELPER
	// this is a recursive function that gets the provenance trail out of the database for a particular provenance id 
	void provenance_recursive_helper( 
		std::vector< std::pair< long long, std::string > >& provenance_trail, ProvenanceID prov_id );
	
private:
	// Session current using
	SessionHandle current_session_;
	
	// Where to save the project
	boost::filesystem::path project_path_;

	// Where the data is being managed
	DataManagerHandle data_manager_;
	
	// Whether the project has changed
	bool changed_;
	
	// The provenance database
	sqlite3* provenance_database_;
	
	bool database_initialized_;
	
};

} // end namespace Seg3D

#endif
