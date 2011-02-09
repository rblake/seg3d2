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

// Core includes
#include <Core/Utils/StringUtil.h>

// Application includes
#include <Application/DatabaseManager/DatabaseManager.h>

namespace Seg3D
{

class DatabaseManagerPrivate {
public:
	sqlite3* database_;
	boost::filesystem::path database_path_;
	std::string error_;
	
};


DatabaseManager::DatabaseManager() :
	private_( new DatabaseManagerPrivate )
{	
	this->private_->error_ = "none";
}

DatabaseManager::~DatabaseManager()
{	
	this->close_database();
}

void DatabaseManager::close_database()
{
	this->load_or_save_database( true );
}

bool DatabaseManager::initialize_database( const boost::filesystem::path& database_path, 
	const std::vector< std::string > database_create_tables_statements )
{
	this->private_->database_path_ = database_path;
	// we check to see if the folder that has been passed as the location for the db actually exists
	if( !boost::filesystem::exists( this->private_->database_path_.parent_path() ) )
	{
		this->private_->error_= "Database could not be created.  The path provided '"
			+ this->private_->database_path_.string() + "', does not exist!";
		return false;
	}
	
	int result;
	
	// we check to see if the database already exists.
	// if it doesn't, we create it using the passed create statements
	if( !boost::filesystem::exists( this->private_->database_path_ ) )
	{
		result = sqlite3_open( ":memory:", &this->private_->database_ );
		if ( result == SQLITE_OK )
		{
			const char* tail;
			sqlite3_stmt* statement;
			
			for( size_t i = 0; i < database_create_tables_statements.size(); ++i )
			{
				statement = NULL;
				sqlite3_prepare_v2( this->private_->database_, database_create_tables_statements[ i ].c_str(), 
				static_cast< int >( database_create_tables_statements[ i ].size() ), &statement, &tail );
				result = sqlite3_step( statement );
				sqlite3_finalize( statement );				

				if(  result != SQLITE_DONE  )
				{
					this->private_->error_=  "Database could not be created.  The create statement provided '"
						+ database_create_tables_statements[ i ] + "', returned error: "
						+ Core::ExportToString( result );
					return false;
				}
			}
		}
		else
		{
			this->private_->error_= "Database could not be created. The database could not be created.";
			return false;
		}
	}
	else
	{
		result = sqlite3_open( ":memory:", &this->private_->database_ );
		if( ( result != SQLITE_OK ) || ( !this->load_or_save_database( false ) ) )
		{
			this->private_->error_=  "Database could not be opened at '"
				+ this->private_->database_path_.string() + "'. And returned error: "
				+ Core::ExportToString( result );
			return false;
		}
	}
	
	return true;
}

bool DatabaseManager::database_query( const std::string& sql_query, ResultSet& results )
{
	int result;
 	const char* tail;
	sqlite3_stmt* statement = NULL;
	sqlite3_prepare_v2( this->private_->database_, sql_query.c_str(), 
		static_cast< int >( sql_query.size() ), &statement, &tail );
	bool at_least_one = false;
	result = sqlite3_step( statement );	
		
	while ( result == SQLITE_ROW )
	{
		at_least_one = true;
	
		std::map< std::string, boost::any > temp_map;
		for( int j = 0; j < sqlite3_column_count( statement ); ++j )
		{
			boost::any temp_any;
			switch( sqlite3_column_type( statement, j ) )
			{
				case SQLITE_TEXT:
				{
					std::string string_result = std::string( (char*)sqlite3_column_text( statement, j ) );
					temp_any = boost::any( string_result );
					break;
				}
				case SQLITE_INTEGER:
				{
					int int_result = sqlite3_column_int( statement, j );
					temp_any = boost::any( int_result );
					break;
				}
			}
			std::string name_result = std::string( sqlite3_column_name( statement, j ) );
			temp_map[ name_result ] = temp_any;
		}
		results.push_back( temp_map );
		result = sqlite3_step( statement );
	}
	
	sqlite3_finalize( statement );
	
	if( !at_least_one && ( result != SQLITE_DONE ) )
	{
		this->private_->error_=  "The select statement had no results. '" + sql_query + "', returned error: "
			+ Core::ExportToString( result );
		return false;
	} 
	
	return true;
}

bool DatabaseManager::database_query_no_return( const std::string& sql_query )
{
	int result;
	const char* tail;
	sqlite3_stmt* statement = NULL;
	sqlite3_prepare_v2( this->private_->database_, sql_query.c_str(), 
		static_cast< int >( sql_query.size() ), &statement, &tail );
	result = sqlite3_step( statement );
	sqlite3_finalize( statement );
	
	if( result != SQLITE_DONE )
	{
		this->private_->error_=  "The insert statement could not be completed.  The insert statement provided '"
			+ sql_query + "', returned error: " + Core::ExportToString( result );
		return false;
	}
	
	return true;
	
}

bool DatabaseManager::database_checkpoint()
{
	return this->load_or_save_database( true );
}

std::string DatabaseManager::get_error()
{
	std::string temp_error = this->private_->error_;
	this->private_->error_ = "none";
	return temp_error;
}

bool DatabaseManager::load_or_save_database( bool is_save )
{
	int result;
	sqlite3* temp_open_database;
	sqlite3* to_database;
	sqlite3* from_database;
	sqlite3_backup* backup_database_object;
	
	result = sqlite3_open( this->private_->database_path_.string().c_str(), &temp_open_database );
	
	if( result != SQLITE_OK ) 
	{
		( void )sqlite3_close( temp_open_database );
		return false;
	}
	
	from_database = ( is_save ? this->private_->database_ : temp_open_database );
	to_database = ( is_save ? temp_open_database : this->private_->database_ );
	
	backup_database_object = sqlite3_backup_init( to_database, "main", from_database, "main" );
	if( backup_database_object )
	{
		( void )sqlite3_backup_step( backup_database_object, -1 );
		( void )sqlite3_backup_finish( backup_database_object );
	}
	
	result = sqlite3_errcode( to_database );
	
	( void )sqlite3_close( temp_open_database );

	return true;	
}



} // end namespace seg3D
