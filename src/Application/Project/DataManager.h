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

#ifndef APPLICATION_DATAMANAGER_DATAMANAGER_H
#define APPLICATION_DATAMANAGER_DATAMANAGER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// STL includes
#include <time.h>

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>


// Core includes
#include <Core/Utils/StringUtil.h>
#include <Core/Utils/Singleton.h>
#include <Core/Utils/Log.h>
#include <Core/Utils/Exception.h>

// Application includes
#include <Core/State/StateHandler.h>


namespace Seg3D
{

// Forward declaration
class DataManager;
	
typedef boost::shared_ptr< DataManager > DataManagerHandle;

// Class definition
class DataManager : public Core::StateHandler
{
	// -- Constructor/Destructor --
public:
	DataManager();
	virtual ~DataManager();
	
protected:
	Core::StateStringVectorHandle sessions_and_datafiles_state_;
	
	Core::StateStringVectorHandle source_datafiles_state_;

public:
	// INITIALIZE:
	// This function loads the values for DataManager from file
	void initialize( const boost::filesystem::path& project_path );
	
	// SAVE_DATAMANAGER_STATE:
	// This function saves the DataManager's state to file
	void save_datamanager_state( const boost::filesystem::path& project_path, 
		const std::string& session_name );

	// REMOVE_SESSION:
	// This function runs when the project emits a signal saying that a session has been deleted,
	// and removes it from the session datafile list
	void remove_session( const std::string& session_name );

	// INSERT_DATA_SOURCE
//	void insert_data_source();

	// GET_SESSION_FILES_VECTOR:
	// this function, when passed the name of a session, returns by reference, a vector of strings
	// containing the list of files that that session depends on.  It returns true if the session
	// exists, and false if the session does not.
	bool get_session_files_vector( const std::string& session_name, 
		std::vector< std::string >& files );

	// CLEAR_DATA_FILE_LIST:
	// funtion clear out the data file list in preperation for starting a new project
	void clear_data_file_list();

	// GET_FILE_SIZE:
	// function that returns the file size calculated from the data files
	long long get_file_size();
	
private:
	// PREP_FOR_SAVE:
	// this function prepares the state variable "essions_and_datafiles_state_" for saving to file
	void prep_for_save( const boost::filesystem::path& project_path, 
		const std::string& session_name );

	// GET_MUTEX:
	// we need a state engine mutex for when we retrieve the list of layers and list of currently
	// used datafiles.
	typedef Core::StateEngine::mutex_type mutex_type;
	typedef Core::StateEngine::lock_type lock_type;
	mutex_type& get_mutex();

private:
	long long disk_space_used;

};

} // end namespace seg3d

#endif // DataManager
