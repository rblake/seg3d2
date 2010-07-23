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

#ifndef CORE_ACTION_ACTIONINFO_H
#define CORE_ACTION_ACTIONINFO_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// STL includes
#include <string>

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/thread.hpp>

namespace Core
{

class ActionInfo;
class ActionInfoPrivate;

typedef boost::shared_ptr< ActionInfo > ActionInfoHandle;
typedef boost::shared_ptr< ActionInfoPrivate > ActionInfoPrivateHandle;

// CLASS ACTIONINFO:
// Information class about an action

class ActionInfo : public boost::noncopyable
{
	// -- constructor --
public:
	ActionInfo( const std::string& definition );
	
public:
	// GET_DEFINITION:
	// Get the definition of the action (in xml format)
	std::string get_definition() const;

	// GET_TYPE:
	// The type/name of the action
	std::string get_type() const;

	// GET_USAGE:
	// The usage string of the action
	std::string get_usage() const;
	
	// GET_DESCRIPTION:
	// Get the description of the action
	std::string get_description() const;
	
	// GET_NUM_ARGUMENTS:
	// The number of arguments for this action
	size_t get_num_arguments() const;
	
	// GET_NUM_KEY_VALUE_PAIRS:
	// THe number of key value pairs associated with this action
	size_t get_num_key_value_pairs() const;
	
	// GET_ARGUMENT:
	// The name of argument <index>
	std::string get_argument( size_t index ) const;

	// GET_ARGUMENT_DESCRIPTION:
	// The description of argument <index>
	std::string get_argument_description( size_t index ) const;
	
	// GET_KEY:
	// The name of key <index>
	std::string get_key( size_t index ) const;
	
	// GET_INDEX:
	// The index of a key
	int get_key_index( const std::string& key ) const;
	
	// GET_DEFAULT_KEY_VALUE:
	// The default value of the key <index>
	std::string get_default_key_value( size_t index ) const;
	
	// GET_KEY_DESCRIPTION:
	// The description of key <index>
	std::string get_key_description( size_t index ) const;
	
private:
	// Private implementation of this class
	ActionInfoPrivateHandle private_;

public:
	typedef boost::mutex mutex_type;
	typedef boost::unique_lock<mutex_type> lock_type;
	
	// GET_MUTEX:
	// Get the mutex that protects the generation of ActionInfo instances
	static mutex_type& GetMutex();

private:
	// the actual mutex that is shared between all the ActionInfo classes
	static mutex_type mutex_;
};

} // end namespace Core

#endif
