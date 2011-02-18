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

// STL includes
#include <time.h>

// Boost includes
#include <boost/filesystem.hpp>

// Core includes
#include <Core/Application/Application.h>
#include <Core/State/StateIO.h>
#include <Core/Utils/StringUtil.h>
#include <Core/Utils/Log.h>
#include <Core/DataBlock/DataBlockManager.h>


// Application includes
#include <Application/Project/Project.h>
#include <Application/PreferencesManager/PreferencesManager.h>
#include <Application/LayerManager/LayerManager.h>

namespace Seg3D
{

Project::Project( const std::string& project_name ) :
	StateHandler( "project", false ),
	changed_( false ),
	database_initialized_( false )
{	
	// Name of the project.
	this->add_state( "project_name", this->project_name_state_, project_name );

	// Whether color from the preferences or from the project are used.
	this->add_state( "save_custom_colors", this->save_custom_colors_state_, false );
	
	// List of sessions stored in this project
	// *** THIS IS DEPRECATED *** // this is only here for backwards compatibility and is not
	// currently used aside from that.
	std::vector< std::string> empty_vector;
	this->add_state( "sessions", this->sessions_state_, empty_vector );
	
	// Notes for the project.
	this->add_state( "project_notes", this->project_notes_state_, empty_vector );

	// Running count of how much the data files consume. It is not exact it is an approximation.
	this->add_state( "project_file_size", this->project_file_size_state_, 0 );
	
	// Name of the session that is currently in use.
	this->add_state( "current_session_name", this->current_session_name_state_, "UnnamedSession" );
	
	// This keeps track of the count for the sessions. Currently the maximum number of sessions 
	// that we handle is 999. 
	this->add_state( "session_count", this->session_count_state_, 0 );
		
	// State of all the 12 colors in the system.	
	this->color_states_.resize( 12 );
	for ( size_t j = 0; j < 12; j++ )
	{
		// Initialize the colors with the default colors from the preference manager
		std::string stateid = std::string( "color_" ) + Core::ExportToString( j );
		this->add_state( stateid, this->color_states_[ j ], 
			PreferencesManager::Instance()->get_default_colors()[ j ] );
	}

	this->add_state( "generation_count", this->generation_count_state_, -1 );

	std::string user;
	Core::Application::Instance()->get_user_name( user );
	
	this->current_session_ = SessionHandle( new Session( "default_session", user ) );
	this->data_manager_ = DataManagerHandle( new DataManager() );

	// Each time an action is execute, check whether it changes the project data, if so
	// mark this in the project, so the UI can query the user for a save action if the application
	// is closed while there is unsaved data.
	this->add_connection( Core::ActionDispatcher::Instance()->post_action_signal_.connect( 
		boost::bind( &Project::set_project_changed, this, _1, _2 ) ) );
}
	
Project::~Project()
{
	this->disconnect_all();
}

bool Project::initialize_from_file( const std::string& project_name )
{
	Core::StateIO stateio;
	if( stateio.import_from_file( this->project_path_ / ( project_name + ".s3d" ) ) &&
		this->load_states( stateio ) )	
	{
		// now lets try and setup the provenance database
		if( !this->create_database_schema() )
		{
			CORE_LOG_ERROR( "Unable to create provenance database!" );
			return false;
		}

		if( !this->sessions_state_->get().empty() )
		{
			this->session_count_state_->set( 0 );
			this->import_old_session_info_into_database();
		}

		// Once we have loaded the state of the sessions_list, we need to validate that the files exist
		this->cleanup_session_database();

		std::string most_recent_session_name;
		if( this->get_most_recent_session_name( most_recent_session_name ) && 
			this->load_session( most_recent_session_name ) )
		{
			this->project_file_size_state_->set( this->get_data_file_size() );
			this->sessions_changed_signal_();
			return true;
		}
	}
	return false;
}

bool Project::load_session( const std::string& session_name )
{
	return  this->current_session_->load( this->project_path_,	session_name );
}

bool Project::delete_session( const std::string& session_name )
{
	boost::filesystem::path session_path = 
		this->project_path_ / "sessions" / ( session_name + ".xml" );

	try 
	{
		boost::filesystem::remove_all( session_path );
	}
	catch(  std::exception& e ) 
	{
		CORE_LOG_WARNING( e.what() );
		return false;
	}

	this->delete_session_from_database( session_name );
	this->delete_unused_data_files();
	this->project_file_size_state_->set( this->get_data_file_size() );

	Core::StateIO stateio;
	stateio.initialize();
	this->save_states( stateio );
	stateio.export_to_file( this->project_path_ / ( this->project_name_state_->get() + ".s3d" ) );

	this->sessions_changed_signal_();
	return true;
}	

bool Project::save_session( const std::string& timestamp, const std::string& session_name )
{
	std::string session_count = Core::ExportToString( this->session_count_state_->get() );
	while( session_count.size() < 3 )
	{
		// we are going to pad session counts so that the os sorts them properly
		session_count = "0" + session_count;
	}

	if( !this->insert_session_into_database( timestamp, session_count + "-" + session_name ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}

	this->current_session_->session_name_state_->set( session_count + "-" + session_name );

	if( !this->current_session_->save( this->project_path_, session_count + "-" + session_name ) )
	{
		return false;
	}

	this->project_file_size_state_->set( this->get_data_file_size() );
	this->sessions_changed_signal_();

	return true;
}

bool Project::save_as( boost::filesystem::path path, const std::string& project_name )
{
	boost::filesystem::path data_path = this->project_path_ / "data";
	boost::filesystem::directory_iterator data_dir_end;
	for( boost::filesystem::directory_iterator data_dir_itr( data_path ); 
		data_dir_itr != data_dir_end; ++data_dir_itr )
	{
		try
		{
			boost::filesystem::copy_file( ( data_path / data_dir_itr->filename() ),
				( path / project_name / "data" / data_dir_itr->filename() ) );
		}
		catch ( std::exception& e ) // any errors that we might get thrown
		{
			CORE_LOG_ERROR( e.what() );
			return false;
		}
	}

	boost::filesystem::path session_path = this->project_path_ / "sessions";
	boost::filesystem::directory_iterator session_dir_end;
	for( boost::filesystem::directory_iterator session_dir_itr( session_path ); 
		session_dir_itr != session_dir_end; ++session_dir_itr )
	{
		try
		{
			boost::filesystem::copy_file( ( session_path / session_dir_itr->filename() ),
				( path / project_name / "sessions"/ session_dir_itr->filename() ) );
		}
		catch ( std::exception& e ) // any errors that we might get thrown
		{
			CORE_LOG_ERROR( e.what() );
			return false;
		}
	}

	return true;
}

bool Project::project_export( boost::filesystem::path path, const std::string& project_name, 
							 const std::string& session_name )
{
	boost::filesystem::path data_path = path / project_name / "data";
	
	boost::filesystem::directory_iterator data_dir_end;
	for( boost::filesystem::directory_iterator data_dir_itr( data_path ); 
		data_dir_itr != data_dir_end; ++data_dir_itr )
	{
		try
		{
			boost::filesystem::copy_file( ( project_path_ / "data" / data_dir_itr->filename() ),
				( data_path / data_dir_itr->filename() ) );
		}
		catch ( std::exception& e ) // any errors that we might get thrown
		{
			CORE_LOG_ERROR( e.what() );
			return false;
		}
	}
	
	try
	{
		boost::filesystem::copy_file( ( project_path_ / "sessions" / ( session_name + ".xml" ) ),
			( path / project_name / "sessions" / ( session_name + ".xml" ) ) );
	}
	catch ( std::exception& e ) // any errors that we might get thrown
	{
		CORE_LOG_ERROR( e.what() );
		return false;
	}

	return true;

}

void Project::set_signal_block( bool on_off )
{
	this->enable_signals( on_off );
}

void Project::set_project_path( const boost::filesystem::path& project_path )
{
	this->project_path_ = project_path;
}

void Project::set_project_changed( Core::ActionHandle action, Core::ActionResultHandle result )
{
	// NOTE: This is executed on the application thread, hence we do not need a lock to read
	// the variable that is only changed on the same thread
	
	if ( this->changed_ == false  && action->changes_project_data() )
	{
		// NOTE: Changing the variable
		Core::Application::lock_type lock( Core::Application::GetMutex() );
		this->changed_ = true;
		return;
	}
}

void Project::reset_project_changed()
{
	Core::Application::lock_type lock( Core::Application::GetMutex() );
	this->changed_ = false;
}

bool Project::check_project_changed()
{
	Core::Application::lock_type lock( Core::Application::GetMutex() );
	return this->changed_;
}

bool Project::pre_save_states( Core::StateIO& state_io )
{
	if( this->save_custom_colors_state_->get() )
	{
		for ( size_t j = 0; j < 12; j++ )
		{
			this->color_states_[ j ]->set( 
				PreferencesManager::Instance()->color_states_[ j ]->get() );
		}
	}
	
	this->generation_count_state_->set( Core::DataBlockManager::Instance()->get_generation_count() );

	return true;
}

bool Project::post_load_states( const Core::StateIO& state_io )
{
	// If the user has chosen to save their custom colors as part of the project, we load them into
	// the preferences manager from the project's state variables
	if( this->save_custom_colors_state_->get() )
	{
		for ( size_t j = 0; j < 12; j++ )
		{
			PreferencesManager::Instance()->color_states_[ j ]->set( 
				this->color_states_[ j ]->get() );
		}
	}

	this->data_manager_->initialize( this->project_path_ );

	// NOTE: Because an earlier version mistakenly did not save this number, it may be initialized
	// to -1, in that case we need to derive it from the actual generation numbers stored in
	// the data directory.
	Core::DataBlock::generation_type generation = this->generation_count_state_->get();
	
	/////////////////////////////////////////////////////////////////////////
	// NOTE: This correction is for backwards compatibility
	//if ( generation == -1 )
	{
		boost::filesystem::path data_path = this->project_path_ / "data";
		if( boost::filesystem::exists( data_path ) )
		{
			boost::filesystem::directory_iterator dir_end;
			for( boost::filesystem::directory_iterator dir_itr( data_path ); 
				dir_itr != dir_end; ++dir_itr )
			{
				if ( boost::filesystem::extension( *dir_itr ) == ".nrrd" )
				{
					Core::DataBlock::generation_type file_generation = -1;
					if ( Core::ImportFromString( boost::filesystem::basename( dir_itr->filename() ),
						file_generation ) )
					{
						if ( file_generation > generation ) generation = file_generation;
					}
				}
			}
		}

		// Ensure the counter will have at least 0
		if ( generation == -1 ) generation = 0;
		generation++;
	}
	/////////////////////////////////////////////////////////////////////////
	
	Core::DataBlockManager::Instance()->set_generation_count( generation );
	this->generation_count_state_->set( generation );

	return true;
}

bool Project::validate_session_name( std::string& session_name )
{
	SessionInfo session;
	if( this->get_session( session, session_name ) )
	{
		if( session.session_name_ != "" )
		{
			return true;
		}
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Provenance Database Functionality //////////////////////////////////
bool Project::create_database_schema()
{
	if( !boost::filesystem::exists( ( this->project_path_ / "provenance" ) ) ) 
	{
		try // to create a project provenance folder because it doesnt exist
		{
			boost::filesystem::create_directory( this->project_path_ / "provenance");
		}
		catch ( std::exception& e ) // any errors that we might get thrown
		{
			CORE_LOG_ERROR( e.what() );
			return false;
		}
	}
	boost::filesystem::path  database_path = this->project_path_ / "provenance" / "project.sqlite";

		std::vector< std::string > create_table_commands;
		create_table_commands.push_back(
			// this table will store the actual provenance step
			"CREATE TABLE provenance_step "
			"(provenance_id INTEGER NOT NULL, "
			"action TEXT NOT NULL, "
			"timestamp TEXT NOT NULL, "
			"user TEXT NOT NULL, "
			"output_number INTEGER NOT NULL, "
			"UNIQUE (provenance_id));" );
			
		create_table_commands.push_back(
			// this is a table relating the input provenance_ids to the output
			"CREATE TABLE provenance_input_relations "
			"(id INTEGER NOT NULL, "
			"provenance_id INTEGER NOT NULL, "
			"input_provenance_id INTEGER NOT NULL, "
			"PRIMARY KEY (id));" );
			
		create_table_commands.push_back(
			// this is a table relating the input provenance_ids to the output
			"CREATE TABLE provenance_output_relations "
			"(id INTEGER NOT NULL, "
			"provenance_id INTEGER NOT NULL, "
			"output_provenance_id INTEGER NOT NULL, "
			"PRIMARY KEY (id));" );
			
		create_table_commands.push_back(
			// this is a table relating the input provenance_ids to the output
			"CREATE TABLE provenance_deleted_relations "
			"(id INTEGER NOT NULL, "
			"provenance_id INTEGER NOT NULL, "
			"deleted_provenance_id INTEGER NOT NULL, "
			"PRIMARY KEY (id));" );
			
		create_table_commands.push_back(
			// this table represents the actions that cause new provenance_id's
			"CREATE TABLE actions "
			"(id INTEGER NOT NULL, "
			"action_id TEXT NOT NULL, "
			"number_of_inputs INTEGER NOT NULL, "
			"number_of_outputs INTEGER NOT NULL, "
			"PRIMARY KEY (id), "
			"UNIQUE (action_id));" );
			
		create_table_commands.push_back(
			// this is a table relating settings to actions
			"CREATE TABLE action_settings "
			"(setting_id INTEGER NOT NULL, "
			"action_id TEXT NOT NULL, "
			"setting_name TEXT NOT NULL, "
			"setting_value TEXT NOT NULL, "
			"PRIMARY KEY (setting_id));" );
			
		create_table_commands.push_back(
			// this table represents the actions that cause new provenance_id's
			"CREATE TABLE sessions "
			"(session_id INTEGER NOT NULL, "
			"session_name TEXT NOT NULL, "
			"username TEXT NOT NULL, "
			"timestamp TEXT NOT NULL, "
			"UNIQUE (session_name),"
			"UNIQUE (timestamp),"
			"PRIMARY KEY (session_id));" );

		create_table_commands.push_back(
			// this is a table relating the input provenance_ids to the output
			"CREATE TABLE data_relations "
			"(id INTEGER NOT NULL, "
			"session_name TEXT NOT NULL, "
			"data_file TEXT NOT NULL, "
			"UNIQUE(session_name, data_file),"
			"PRIMARY KEY (id));" );
		

	if(	this->initialize_database( database_path, create_table_commands ) )
	{
		this->database_initialized_ = true;
		return true;
	}
	
	return false;
}

// This function is mostly just a placeholder.  Currently it just registers the actions.  We will probably want to 
// create a Providence Object and then add it to the db.
bool Project::add_to_provenance_database( ProvenanceStepHandle& step )
{
	// Print diagnostics
	//step->print();
	
	std::string action_desc = step->get_action();
	std::string user_name = step->get_username();
	std::string timestamp = step->get_timestamp();
		
	ProvenanceIDList output_list = step->get_output_provenance_ids();
	ProvenanceIDList input_list = step->get_input_provenance_ids();
	ProvenanceIDList deleted_list = step->get_deleted_provenance_ids();
	
	std::string insert_statement;
	
	for( size_t i = 0; i < output_list.size(); ++i )
	{
		insert_statement = 
			"INSERT INTO provenance_step (provenance_id, action, timestamp, user, output_number) VALUES(" + 
			Core::ExportToString( output_list[ i ] )+ ", '" + action_desc + "', '" + 
			timestamp + "', '" + user_name + "', " + Core::ExportToString( i+1 )+ ");";
		
		if( !this->database_query_no_return( insert_statement ) )
		{
			CORE_LOG_ERROR( this->get_error() );
			return false;
		}
		
		for( size_t j = 0; j < input_list.size(); ++j )
		{
			insert_statement = "INSERT INTO provenance_input_relations (provenance_id, input_provenance_id) VALUES(" + 
				Core::ExportToString( output_list[ i ] ) + ", " +
				Core::ExportToString( input_list[ j ] ) + ");";

			if( !this->database_query_no_return( insert_statement ) )
			{
				CORE_LOG_ERROR( this->get_error() );
				return false;
			}
		}

		for( size_t j = 0; j < output_list.size(); ++j )
		{
			insert_statement = "INSERT INTO provenance_output_relations (provenance_id, output_provenance_id) VALUES(" + 
				Core::ExportToString( output_list[ i ] ) + ", " +
				Core::ExportToString( output_list[ j ] ) + ");";

			if( !this->database_query_no_return( insert_statement ) )
			{
				CORE_LOG_ERROR( this->get_error() );
				return false;
			}
		}

		for( size_t j = 0; j < deleted_list.size(); ++j )
		{
			insert_statement = "INSERT INTO provenance_deleted_relations "
				"(provenance_id, deleted_provenance_id) VALUES(" + 
				Core::ExportToString( output_list[ i ] ) + ", " +
				Core::ExportToString( deleted_list[ j ] ) + ");";

			this->database_query_no_return( insert_statement );
		}	
		
	}
	return true;	
}

void Project::close_provenance_database()
{
	this->close_database();
}

void Project::checkpoint_provenance_database()
{
	if( !this->database_initialized_ ) return;
	this->database_checkpoint();
}

bool Project::insert_session_into_database( const std::string& timestamp, 
	const std::string& session_name, const std::string& user_name /*= "" */ )
{
	std::string user;
	bool importing_old_sessions = false;
	
	// if we have a username it means that we are importing old sessions from file rather than saving
	// the current session 
	if( user_name == "" )
	{
		Core::Application::Instance()->get_user_name( user );
	}
	else
	{
		user = user_name;
		importing_old_sessions = true;
	}
		
	std::string insert_statement = "INSERT OR IGNORE INTO sessions "
		"(session_name, username, timestamp) VALUES('" + 
		session_name+ "', '" + user + "', '" + timestamp + "');";

	if( !this->database_query_no_return( insert_statement ) )
	{
		CORE_LOG_WARNING( "Chill out. Wait at least a second before doing another save..." );
		return false;
	}
	
	if( !importing_old_sessions )
	{
		std::vector< LayerHandle > current_layers;
		LayerManager::Instance()->get_layers( current_layers );

		for( size_t i = 0; i < current_layers.size(); ++i )
		{
			Core::DataBlock::generation_type generation = current_layers[ i ]->get_generation();

			if( generation == -1 ) continue;

			insert_statement = "INSERT OR IGNORE INTO data_relations "
				"(session_name, data_file) VALUES('" + session_name + 
				"', '" + Core::ExportToString( generation ) +	".nrrd');";

			if( !this->database_query_no_return( insert_statement ) )
			{
				CORE_LOG_ERROR( this->get_error() );
				return false;
			}
		}
	}
	
	
	this->session_count_state_->set( this->session_count_state_->get() + 1 );
	return true;
}

bool Project::delete_session_from_database( const std::string& session_name )
{
	std::string delete_statement = "DELETE FROM sessions WHERE (session_name = '" + 
		session_name + "');";

	if( !this->database_query_no_return( delete_statement ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}
	
	delete_statement = "DELETE FROM data_relations WHERE (session_name = '" + 
		session_name + "');";

	if( !this->database_query_no_return( delete_statement ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}
	
	return true;
}

bool Project::get_all_sessions( std::vector< SessionInfo >& sessions )
{
	ResultSet result_set;
	std::string select_statement = "SELECT * FROM sessions ORDER BY session_id DESC;";
	if( !this->database_query( select_statement, result_set ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}

	for( size_t i = 0; i < result_set.size(); ++i )
	{
		sessions.push_back( SessionInfo( 
			boost::any_cast< std::string >( ( result_set[ i ] )[ "session_name" ] ),
			boost::any_cast< std::string >( ( result_set[ i ] )[ "username" ] ),
			boost::any_cast< std::string >( ( result_set[ i ] )[ "timestamp" ] ) ) );
	}

	return true;
}

bool Project::get_session( SessionInfo& session, const std::string& session_name )
{
	ResultSet result_set;
	std::string select_statement = "SELECT * FROM sessions WHERE session_name ='" + 
		session_name + "';";
	if( !this->database_query( select_statement, result_set ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}

	session = SessionInfo( 
		boost::any_cast< std::string >( ( result_set[ 0 ] )[ "session_name" ] ),
		boost::any_cast< std::string >( ( result_set[ 0 ] )[ "username" ] ),
		boost::any_cast< std::string >( ( result_set[ 0 ] )[ "timestamp" ] ) );

	return true;
}

bool Project::get_most_recent_session_name( std::string& session_name )
{
	ResultSet result_set;
	std::string select_statement = "SELECT * FROM sessions ORDER BY session_id DESC;";
	if( !this->database_query( select_statement, result_set ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}

	session_name = boost::any_cast< std::string >( ( result_set[ 0 ] )[ "session_name" ] );

	return true;
}

void Project::cleanup_session_database()
{
	std::vector< SessionInfo > sessions;
	if( get_all_sessions( sessions ) )
	{
		for( size_t i = 0; i < sessions.size(); ++i )
		{
			boost::filesystem::path session_path = this->project_path_ / "sessions" 
				/ ( sessions[ i ].session_name_ + ".xml" );
			if( !boost::filesystem::exists( session_path ) )
			{
				this->delete_session_from_database( sessions[ i ].session_name_ );
			}
		}
	}
}

void Project::import_old_session_info_into_database()
{
	std::vector< std::string > session_vector = this->sessions_state_->get();
	int number_of_sessions_to_import = static_cast< int >( session_vector.size() );
	for( int i = ( number_of_sessions_to_import - 1 ); i >= 0; i-- )
	{

		std::string old_session_name;
		Core::ImportFromString( session_vector[ i ], old_session_name );
		
		boost::filesystem::path  old_path = this->project_path_ / "sessions" / 
			( old_session_name + ".xml" );
		
		std::vector< std::string > old_session_name_vector = 
			Core::SplitString( old_session_name, " - " );
		
		std::string session_count = Core::ExportToString( this->session_count_state_->get() );
		while( session_count.size() < 3 )
		{
			// we are going to pad session counts so that the os sorts them properly
			session_count = "0" + session_count;
		}
		
		std::string new_session_name = session_count + "-" + old_session_name_vector[ 1 ];
		
		boost::filesystem::path  new_path = this->project_path_ / "sessions" / 
			( new_session_name + ".xml" );
			
		if( !boost::filesystem::exists( old_path ) ) continue;	
		boost::filesystem::rename( old_path, new_path );
		
		std::string day = ( Core::SplitString( old_session_name_vector[ 0 ], "-" ) )[ 0 ];
		std::string month = ( Core::SplitString( old_session_name_vector[ 0 ], "-" ) )[ 1 ];
		std::string year = ( Core::SplitString( old_session_name_vector[ 0 ], "-" ) )[ 2 ];
		
		std::string hour = ( Core::SplitString( old_session_name_vector[ 0 ], "-" ) )[ 3 ];
		std::string minute = ( Core::SplitString( old_session_name_vector[ 0 ], "-" ) )[ 4 ];
		std::string second = ( Core::SplitString( old_session_name_vector[ 0 ], "-" ) )[ 5 ];
		
		std::string timestamp = year + "-" + month + "-" + day + " " + 
			hour + ":" + minute + ":" + second;
		
		// Here we insert the converted session name into the database
		bool success = this->insert_session_into_database( timestamp, new_session_name, 
			old_session_name_vector[ 2 ] );
		
		// Now if it had data associated with it, put it in the database too.
		std::vector< std::string > files;
		if( this->data_manager_->get_session_files_vector( old_session_name, files ) )
		{
			for( size_t i = 0; i < files.size(); ++i )
			{
				std::string insert_statement = "INSERT OR IGNORE INTO data_relations "
					"(session_name, data_file) VALUES('" + new_session_name + 
					"', '" + files[ i ] + "');";

				if( !this->database_query_no_return( insert_statement ) )
				{
					CORE_LOG_ERROR( this->get_error() );
					return;
				}
			}
		}
		
	}
	
	// Finally we set the session list to an empty vector so we aren't tempted to do this again.
	std::vector< std::string > empty_vector;
	this->sessions_state_->set( empty_vector );
}

long long Project::get_data_file_size()
{
	ResultSet result_set;
	std::string select_statement = "SELECT DISTINCT data_file FROM data_relations;";
	if( !this->database_query( select_statement, result_set ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return false;
	}

	long long size = 0;
	
	for( size_t i = 0; i < result_set.size(); ++i )
	{
		boost::filesystem::path data_file_path = this->project_path_ / "data" / 
			boost::any_cast< std::string >( ( result_set[ i ] )[ "data_file" ] );
		
		if( boost::filesystem::exists( data_file_path ) )
		{
			size += static_cast< long long >( boost::filesystem::file_size( data_file_path ) );
		}
	}
	
	return size;
}

void Project::delete_unused_data_files()
{
	ResultSet result_set;
	std::string select_statement = "SELECT DISTINCT data_file FROM data_relations;";
	if( !this->database_query( select_statement, result_set ) )
	{
		CORE_LOG_ERROR( this->get_error() );
		return;
	}
	
	boost::filesystem::path path =  this->project_path_ / "data";
	
	std::vector< std::string > existing_data_files;
	
	boost::filesystem::directory_iterator dir_end;
	for( boost::filesystem::directory_iterator dir_itr( path ); 
		dir_itr != dir_end; ++dir_itr )
	{
		existing_data_files.push_back( dir_itr->filename() );
	}
	
	for( size_t j = 0; j < existing_data_files.size(); ++j )
	{
		bool found = false;
		for( size_t i = 0; i < result_set.size(); ++i )
		{
			if( boost::any_cast< std::string >( ( result_set[ i ] )[ "data_file" ] ) == 
				existing_data_files[ j ] ) 
			{
				found = true;
				break;	
			}
		}
		if( !found )
		{
			boost::filesystem::remove( path / existing_data_files[ j ] );
		}
	}	
}



} // end namespace Seg3D
