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

#ifndef APPLICATION_DATABASEMANAGER_DATABASEMANAGER_H
#define APPLICATION_DATABASEMANAGER_DATABASEMANAGER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

#include <map>

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>

// Sqlite includes
#include <Externals/sqlite/sqlite3.h>


namespace Seg3D
{

typedef std::vector< std::map< std::string, boost::any > > ResultSet;

// Forward declaration
class DatabaseManagerPrivate;

// Class definition
class DatabaseManager 
{
	// -- Constructor/Destructor --
public:
	DatabaseManager();
	virtual ~DatabaseManager();
	
public:
	// RUN_SQL_STATEMENT:
	// Execute the given SQL statement on the database. If the statement generates
	// any results, they will be put in the result set.
	// Returns true on success, otherwise false.
	bool run_sql_statement( const std::string& sql_str, ResultSet& results );

	// RUN_SQL_STATEMENT:
	// Execute the given SQL statement on the database.
	// Returns true on success, otherwise false.
	bool run_sql_statement( const std::string& sql_str );
	
	// DATABASE_CHECKPOINT:
	// this is for writing the database to disk.  We need to do this because the db is stored in
	// memory	
	bool database_checkpoint();
	
	// GET_ERROR:
	// this is an acessor function for getting the error if one of the query's returned false
	std::string get_error();
	
	// CLOSE_DATABASE:
	// this function writes the database to disk and closes the database 
	void close_database();
	
	// CREATE_DATABASE_SCHEMA:
	// this function call the intitialize database function
	virtual bool create_database_schema(){ return false; }
	
	// INITIALIZE_DATABASE:
	// this function does the actual creation of the databases 
	bool initialize_database( const boost::filesystem::path& database_path, 
		const std::vector< std::string > database_create_tables_statements );
	
private:
	// LOAD_OR_SAVE_DATABASE:
	// function that handles both loading or saving of the db to file
	bool load_or_save_database( bool is_save );

private:
	boost::shared_ptr< DatabaseManagerPrivate > private_;

};

} // end namespace seg3d

#endif

