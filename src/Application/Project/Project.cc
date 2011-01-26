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

// Application includes
#include <Application/Project/Project.h>
#include <Application/PreferencesManager/PreferencesManager.h>

namespace Seg3D
{

Project::Project( const std::string& project_name ) :
	StateHandler( "project", false ),
	changed_( false )
{	
	this->add_state( "project_name", this->project_name_state_, project_name );
	this->add_state( "save_custom_colors", this->save_custom_colors_state_, false );
	
	std::vector< std::string> empty_vector;
	this->add_state( "sessions", this->sessions_state_, empty_vector );
	this->add_state( "project_notes", this->project_notes_state_, empty_vector );
	this->add_state( "project_file_size", this->project_file_size_state_, 0 );
	this->add_state( "current_session_name", this->current_session_name_state_, "UnnamedSession" );
		
	this->color_states_.resize( 12 );
	for ( size_t j = 0; j < 12; j++ )
	{
		std::string stateid = std::string( "color_" ) + Core::ExportToString( j );
		this->add_state( stateid, this->color_states_[ j ], 
			PreferencesManager::Instance()->get_default_colors()[ j ] );
	}

	this->current_session_ = SessionHandle( new Session( "default_session" ) );
	this->data_manager_ = DataManagerHandle( new DataManager() );

	this->add_connection( Core::ActionDispatcher::Instance()->post_action_signal_.connect( 
		boost::bind( &Project::set_project_changed, this, _1, _2 ) ) );
		
	this->add_connection( Core::ActionDispatcher::Instance()->post_action_signal_.connect( 
		boost::bind( &Project::add_to_provenance_database, this, _1, _2 ) ) );
}
	
Project::~Project()
{
	this->disconnect_all();
}



void Project::set_project_changed( Core::ActionHandle action, Core::ActionResultHandle result )
{
	// NOTE: This is executed on the application thread, hence we do not need a lock to read
	// the variable that is only changed on the same thread
	
	if ( action->changes_project_data() && this->changed_ == false )
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

bool Project::initialize_from_file( const std::string& project_name )
{
	Core::StateIO stateio;
	if( stateio.import_from_file( this->project_path_ / ( project_name + ".s3d" ) ) &&
		this->load_states( stateio ) )	
	{
		if( this->load_session( this->get_session_name( 0 ) ) )
		{
			return true;
		}
	}
	return false;
}
	
bool Project::load_session( const std::string& session_name )
{
	return  this->current_session_->load( this->project_path_,	session_name );
}

bool Project::save_session( const std::string& session_name )
{
	this->current_session_->session_name_state_->set( session_name );

	if( this->current_session_->save( this->project_path_, session_name ) )
	{
		this->add_session_to_list( session_name );
		return true;
	}
	return false;
}

bool Project::delete_session( const std::string& session_name )
{
	boost::filesystem::path session_path = 
		this->project_path_ / "sessions" / ( session_name + ".xml" );
	
	std::vector< std::string > temp_sessions_vector = this->sessions_state_->get();
	
	for( int i = 0; i < static_cast< int >( temp_sessions_vector.size() ); ++i )
	{
		if( temp_sessions_vector[ i ] == session_name )
		{
			temp_sessions_vector.erase( temp_sessions_vector.begin() + i );
			break;
		}
	}
	
	this->sessions_state_->set( temp_sessions_vector );
	
	try 
	{
		boost::filesystem::remove_all( session_path );
	}
	catch(  std::exception& e ) 
	{
		CORE_LOG_WARNING( e.what() );
		return false;
	}
	
	this->data_manager_->remove_session( session_name );
	
	Core::StateIO stateio;
	stateio.initialize();
	this->save_states( stateio );
	stateio.export_to_file( this->project_path_ / ( this->project_name_state_->get() + ".s3d" ) );
	
	return true;
}	
	
void Project::add_session_to_list( const std::string& session_name )
{
	std::vector< std::string > temp_sessions_vector = this->sessions_state_->get();
	
	for( int i = 0; i < static_cast< int >( temp_sessions_vector.size() ); ++i )
	{
		if( ( temp_sessions_vector[ i ] == "" ) || ( temp_sessions_vector[ i ] == "]" ) )
		{
			continue;
		}
		
		if( ( ( Core::SplitString( temp_sessions_vector[ i ], " - " ) )[ 1 ]== "AutoSave" ) && 
			( ( Core::SplitString( session_name, " - " ) )[ 1 ]== "AutoSave" ) )
		{
			this->delete_session( temp_sessions_vector[ i ] );
			temp_sessions_vector.erase( temp_sessions_vector.begin() + i );
			break;
		}
	}
	temp_sessions_vector.insert( temp_sessions_vector.begin(), ( session_name ) );
	this->sessions_state_->set( temp_sessions_vector );

	this->data_manager_->save_datamanager_state( this->project_path_, session_name );

	this->project_file_size_state_->set( this->data_manager_->get_file_size() );
	
}

std::string Project::get_session_name( int index )
{
	std::vector< std::string > temp_sessions_vector = this->sessions_state_->get();
	std::string session_name = "";
	for( int i = 0; i < static_cast< int >( temp_sessions_vector.size() ); ++i )
	{
		if( i == index )
		{
			session_name = temp_sessions_vector[ i ];
			break;
		}
	}
	return session_name;
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
	// Once we have loaded the state of the sessions_list, we need to validate that the files exist
	this->cleanup_session_list();

	return true;
}

bool Project::validate_session_name( std::string& session_name )
{
	std::vector< std::string > temp_sessions_vector = this->sessions_state_->get();
	for( int i = 0; i < static_cast< int >( temp_sessions_vector.size() ); ++i )
	{
		if( temp_sessions_vector[ i ] == ( session_name ) )
		{
			return true;
		}
	}
	return false;
}

void Project::set_project_path( const boost::filesystem::path& project_path )
{
	this->project_path_ = project_path;
}

void Project::cleanup_session_list()
{
	std::vector< std::string > sessions_vector = this->sessions_state_->get();
	std::vector< std::string > new_session_vector;

	for( int i = 0; i < static_cast< int >( sessions_vector.size() ); ++i )
	{
		if( ( sessions_vector[ i ] != "" ) && ( sessions_vector[ i ] != "]" ) )
		{
			boost::filesystem::path session_path = this->project_path_ / "sessions" 
				/ ( sessions_vector[ i ] + ".xml" );
			if( boost::filesystem::exists( session_path ) )
			{
				new_session_vector.push_back( sessions_vector[ i ] );
			}
			else
			{
				this->data_manager_->remove_session( sessions_vector[ i ] );
			}
		}
	}
	this->sessions_state_->set( new_session_vector );
}

bool Project::project_export( boost::filesystem::path path, const std::string& project_name, 
	const std::string& session_name )
{
	this->data_manager_->save_datamanager_state( path, session_name );

	std::vector< std::string > file_vector;
	if( this->data_manager_->get_session_files_vector( session_name, file_vector ) )
	{
		for( int i = 0; i < static_cast< int >( file_vector.size() ); ++i )
		{
			if( !boost::filesystem::exists( path / project_name / "data" / file_vector[ i ] ) )
			{
				try
				{
					boost::filesystem::copy_file( ( project_path_ / "data" / file_vector[ i ] ),
						( path / project_name / "data" / file_vector[ i ] ) );
				}
				catch ( std::exception& e ) // any errors that we might get thrown
				{
					CORE_LOG_ERROR( e.what() );
					return false;
				}
			}
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

////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Provenance Database Functionality //////////////////////////////////

bool Project::create_database_scheme()
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
	boost::filesystem::path  database_path = this->project_path_ / "provenance" / "provenancedatabase.sqlite";

	if( boost::filesystem::exists( database_path ) ) return true;

	lock_type lock( this->get_mutex() );
	if ( sqlite3_open( database_path.string().c_str(), &this->provenance_database_ ) == SQLITE_OK )
	{
		std::vector< std::string > create_table_commands;
		create_table_commands.push_back(
			// this table will store the actual provenance step
			"CREATE TABLE provenance_step "
			"(provenance_id INTEGER NOT NULL, "
			"output_sequence_number INTEGER NOT NULL, "
			"action_id VARCHAR(255) NOT NULL, "
			"user VARCHAR(255) NOT NULL, "
			"mask_bit INTEGER NOT NULL, "
			"generation_number INTEGER NOT NULL, "
			"PRIMARY KEY (provenance_id));" );
			
		create_table_commands.push_back(
			// this is a table relating the input provenance_ids to the output
			"CREATE TABLE provenance_relations "
			"(id INTEGER NOT NULL, "
			"output_provenance_id INTEGER NOT NULL, "
			"input_provenance_id INTEGER NOT NULL, "
			"PRIMARY KEY (id));" );
		
		create_table_commands.push_back(
			// this table represents the actions that cause new provenance_id's
			"CREATE TABLE actions "
			"(id INTEGER NOT NULL, "
			"action_id VARCHAR(255) NOT NULL, "
			"number_of_inputs INTEGER NOT NULL, "
			"number_of_outputs INTEGER NOT NULL, "
			"PRIMARY KEY (id), "
			"UNIQUE (action_id));" );
			
		create_table_commands.push_back(
			// this is a table relating settings to actions
			"CREATE TABLE action_settings "
			"(setting_id INTEGER NOT NULL, "
			"action_id VARCHAR(255) NOT NULL, "
			"setting_name VARCHAR(255) NOT NULL, "
			"setting_value VARCHAR(255) NOT NULL, "
			"PRIMARY KEY (setting_id));" );
		
		for( unsigned int i = 0; i < create_table_commands.size(); ++i )
		{
			char** error = NULL;
			sqlite3_exec( this->provenance_database_, create_table_commands[i].c_str(), NULL, 0, error );
		}
		sqlite3_close( this->provenance_database_ );
	}
	
	// at this point the database file should be created, so we check to make sure they have.
	if( !boost::filesystem::exists( database_path ) ) return false;
	
	return true;
}

// This function is mostly just a placeholder.  Currently it just registers the actions.  We will probably want to 
// create a Providence Object and then add it to the db.
void  Project::add_to_provenance_database( Core::ActionHandle action, Core::ActionResultHandle result )
{
	this->register_action_in_database( action );
	return;

}

// This function simple registers actions in the database as they come in.  It will only put one 
// entry per actiontype.
bool Project::register_action_in_database( Core::ActionHandle action )
{
	lock_type lock( this->get_mutex() );
	bool success = false;
	
	boost::filesystem::path  database_path = this->project_path_ / "provenance" / "provenancedatabase.sqlite";

	if ( sqlite3_open( database_path.string().c_str(), 
		&this->provenance_database_ ) == SQLITE_OK )
	{
		std::string action_name = action->get_type();
		
		const char* tail;
		sqlite3_stmt* statement = NULL;
		std::string insert_statement = "INSERT INTO actions (action_id, number_of_inputs, number_of_outputs) VALUES(?, ?, ?)";
		sqlite3_prepare_v2( this->provenance_database_, insert_statement.c_str(), 
			static_cast< int >( insert_statement.size() ), &statement, &tail );
		sqlite3_bind_text( statement, 1, action_name.c_str(), 
			static_cast< int >( action_name.size() ), SQLITE_TRANSIENT );
		sqlite3_bind_int( statement, 2, 1 );
		sqlite3_bind_int( statement, 3, 1 );
		success = ( sqlite3_step( statement ) == SQLITE_DONE );
		sqlite3_finalize( statement );
		sqlite3_close( this->provenance_database_ );

	}
	
	this->check_database();
	return success;
}

void Project::check_database()
{
	lock_type lock( this->get_mutex() );
	bool success = false;
	int counter = 0;
	boost::filesystem::path  database_path = this->project_path_ / "provenance" / "provenancedatabase.sqlite";

	if ( sqlite3_open(  database_path.string().c_str(), 
		&this->provenance_database_ ) == SQLITE_OK )
	{
		const char* tail;
		sqlite3_stmt* statement = NULL;
		std::string select_statement = "SELECT * FROM actions ORDER BY id DESC";
		sqlite3_prepare_v2( this->provenance_database_, select_statement.c_str(), 
			static_cast< int >( select_statement.size() ), &statement, &tail );
		int i = 0;
		while ( sqlite3_step( statement ) == SQLITE_ROW )
		{
			std::string action_name =std::string( (char*)sqlite3_column_text( statement, 1 ) );
			success = true;
			counter++;
		} 
		sqlite3_finalize( statement );
		sqlite3_close( this->provenance_database_ );
	}
}



} // end namespace Seg3D
