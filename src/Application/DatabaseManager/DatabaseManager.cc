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

// Core includes
#include <Core/Utils/Lockable.h>


namespace Seg3D
{

class DatabaseManagerPrivate : public Core::RecursiveLockable {
public:
	// The actual database
	sqlite3* database_;
	
	// Where the database is located on disk
	boost::filesystem::path database_path_;
};


DatabaseManager::DatabaseManager() :
	private_( new DatabaseManagerPrivate )
{
	// Open a default database
	int result = sqlite3_open( ":memory:", &this->private_->database_ );
	if ( result != SQLITE_OK )
	{
		this->private_->database_ = 0;
	}
}


DatabaseManager::~DatabaseManager()
{	
	// We need to close the database to avoid memory being leaked
	if ( this->private_->database_ )
	{
		sqlite3_close( this->private_->database_ );
	}
}

bool DatabaseManager::create_database( const std::vector< std::string >& create_tables_statements,
	std::string& error )
{
	const char* tail;
	sqlite3_stmt* statement;
	
	for( size_t i = 0; i < create_tables_statements.size(); ++i )
	{
		statement = 0;
		sqlite3_prepare_v2( this->private_->database_, create_tables_statements[ i ].c_str(), 
			static_cast< int >( create_tables_statements[ i ].size() ), &statement, &tail );
		int result = sqlite3_step( statement );
		sqlite3_finalize( statement );				

		if( result != SQLITE_DONE  )
		{
			error =  "Database could not be created. The create statement provided '"
				+ create_tables_statements[ i ] + "', returned error: "
				+ Core::ExportToString( result );
			return false;
		}
	}
	
	return true;
}


bool DatabaseManager::run_sql_statement( const std::string& sql_str, std::string& error )
{
	ResultSet dummy_results;
	return this->run_sql_statement( sql_str, dummy_results, error );
}


bool DatabaseManager::run_sql_statement( const std::string& sql_str, ResultSet& results, 
	std::string& error )
{
	DatabaseManagerPrivate::lock_type lock( this->private_->get_mutex() );
	int result;
	const char* tail;
	sqlite3_stmt* statement = NULL;
	sqlite3_prepare_v2( this->private_->database_, sql_str.c_str(), 
		static_cast< int >( sql_str.size() ), &statement, &tail );

	while ( ( result = sqlite3_step( statement ) ) == SQLITE_ROW )
	{
		std::map< std::string, boost::any > temp_map;
		for( int j = 0; j < sqlite3_column_count( statement ); ++j )
		{
			boost::any temp_any;
			switch( sqlite3_column_type( statement, j ) )
			{
			case SQLITE_TEXT:
			case SQLITE_BLOB:
				{
					std::string string_result = std::string( reinterpret_cast< const char* >( 
						sqlite3_column_text( statement, j ) ) );
					temp_any = boost::any( string_result );
					break;
				}
			case SQLITE_INTEGER:
				{
					long long int_result = sqlite3_column_int64( statement, j );
					temp_any = boost::any( int_result );
					break;
				}
			case SQLITE_FLOAT:
				{
					double double_result = sqlite3_column_double( statement, j );
					temp_any = boost::any( double_result );
					break;
				}
			case SQLITE_NULL:
			default:
				break;
			}
			std::string column_name = std::string( sqlite3_column_name( statement, j ) );
			temp_map[ column_name ] = temp_any;
		}
		results.push_back( temp_map );
	}

	sqlite3_finalize( statement );

	if( result != SQLITE_DONE )
	{
		error =  "The SQL statement '" + sql_str + "' returned error code: "
			+ Core::ExportToString( result );
		return false;
	} 

	return true;
}


bool DatabaseManager::load_database( const boost::filesystem::path& database_file, 
	std::string& error )
{
	DatabaseManagerPrivate::lock_type lock( this->private_->get_mutex() );

	int result;
	sqlite3* temp_open_database;
	sqlite3_backup* backup_database_object;
	
	result = sqlite3_open( database_file.string().c_str(), &temp_open_database );
	
	if ( result != SQLITE_OK ) 
	{
		sqlite3_close( temp_open_database );
		error = std::string( "Could not open database file '" ) + database_file.string() + "'.";
		return false;
	}
	
	backup_database_object = 
		sqlite3_backup_init( this->private_->database_, "main", temp_open_database, "main" );
	
	if ( backup_database_object )
	{
		sqlite3_backup_step( backup_database_object, -1 );
		sqlite3_backup_finish( backup_database_object );
	}
	
	result = sqlite3_errcode( this->private_->database_ );
	if ( result != SQLITE_OK ) 
	{
		error = "Internal error in database.";
		return false;
	}
	
	sqlite3_close( temp_open_database );

	error = "";
	return true;	
}


bool DatabaseManager::save_database( const boost::filesystem::path& database_file, 
	std::string& error )
{
	DatabaseManagerPrivate::lock_type lock( this->private_->get_mutex() );
	int result;
	sqlite3* temp_open_database;
	sqlite3_backup* backup_database_object;
	
	result = sqlite3_open( database_file.string().c_str(), &temp_open_database );
	
	if ( result != SQLITE_OK ) 
	{
		sqlite3_close( temp_open_database );
		error = std::string( "Could not open database file '" ) + database_file.string() + "'.";
		return false;
	}
	
	backup_database_object = 
		sqlite3_backup_init( temp_open_database, "main", this->private_->database_, "main" );
	
	if ( backup_database_object )
	{
		sqlite3_backup_step( backup_database_object, -1 );
		sqlite3_backup_finish( backup_database_object );
	}
	
	result = sqlite3_errcode( temp_open_database );
	if ( result != SQLITE_OK ) 
	{
		error = "Internal error in database.";
		return false;
	}
	
	sqlite3_close( temp_open_database );

	error = "";
	return true;	
}


} // end namespace seg3D
