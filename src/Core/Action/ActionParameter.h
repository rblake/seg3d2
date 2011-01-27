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

#ifndef CORE_ACTION_ACTIONPARAMETER_H
#define CORE_ACTION_ACTIONPARAMETER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// STL
#include <string>
#include <algorithm>

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

// Core
#include <Core/Utils/Log.h>
#include <Core/Utils/StringUtil.h>

namespace Core
{

// PARAMETER CLASSES
// ACTIONPARAMETER<TYPE> and ACTIONPARAMETERVARIANT
// These two classes are both used to store parameters the first one explicitly 
// states the type and the second one does type conversion. The difference is
// that for the first class type the compiler checks the type integrity and 
// hence is the preferred model. However types are not always know up front
// or the action is so general that the parameter is not of a specied class
// for the latter case do we use the variant version.


// ACTIONPARAMETERBASE:
// Base class needed for uniform access to import and export the value
// in a uniform way.

class ActionParameterBase;
typedef boost::shared_ptr< ActionParameterBase > ActionParameterBaseHandle;

class ActionParameterBase
{
	// -- define handle --
public:
	typedef boost::shared_ptr< ActionParameterBase > Handle;

	// -- destructor --
public:
	virtual ~ActionParameterBase();

	// -- functions for accessing data --
	// EXPORT_TO_STRING
	// export the contents of the parameter to string
	virtual std::string export_to_string() const = 0;

	// IMPORT_FROM_STRING
	// import a parameter from a string. The function returns true
	// if the import succeeded
	virtual bool import_from_string( const std::string& str ) = 0;

};

// ACTIONPARAMETER:
// Parameter for an action.

// Forward declaration:
template< class T > class ActionParameter;

// Class definition:
template< class T >
class ActionParameter : public ActionParameterBase
{

	// -- define handle --
public:
	typedef boost::shared_ptr< ActionParameter< T > > Handle;

	// -- constructor/destructor --
public:
	ActionParameter( T* value_ptr ) :
	  value_ptr_( value_ptr )
	{
	}

	virtual ~ActionParameter()
	{
	}

	// -- access to value --
public:

	// EXPORT_TO_STRING
	// export the contents of the parameter to string
	virtual std::string export_to_string() const
	{
		return ExportToString( *( this->value_ptr_ ) );
	}

	// IMPORT_FROM_STRING
	// import a parameter from a string. The function returns true
	// if the import succeeded
	virtual bool import_from_string( const std::string& str )
	{
		return ( ImportFromString( str, *( this->value_ptr_ ) ) );
	}

private:
	// The actual value
	T* value_ptr_;
};


} // namespace Core


#endif
