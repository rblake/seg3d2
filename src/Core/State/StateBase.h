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

#ifndef CORE_STATE_STATEBASE_H
#define CORE_STATE_STATEBASE_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// STL includes
#include <string>

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

// Action includes
#include <Core/Action/ActionContext.h>
#include <Core/Action/ActionParameter.h>

namespace Core
{

// -- STATEBASE --

class StateBase;
typedef boost::shared_ptr< StateBase > StateBaseHandle;
typedef boost::weak_ptr< StateBase > StateBaseWeakHandle;

class StateBase : public boost::noncopyable
{

	// -- constructor / destructor --
public:
	StateBase(const std::string& stateid) :
		stateid_( stateid )
	{
	}
	
	virtual ~StateBase()
	{
	}

	// -- functions for accessing data --

	// EXPORT_TO_STRING:
	// Convert the contents of the State into a string
	virtual std::string export_to_string() const = 0;

	// IMPORT_FROM_STRING:
	// Set the State from a string
	virtual bool import_from_string( const std::string& str, Core::ActionSource source =
		Core::ActionSource::NONE_E ) = 0;

protected:
	friend class ActionSet;
	friend class ActionGet;
	friend class StateHandler;

	// EXPORT_TO_VARIANT
	// Export the state data to a variant parameter
	virtual void export_to_variant( Core::ActionParameterVariant& variant ) const = 0;

	// IMPORT_FROM_VARIANT:
	// Import the state data from a variant parameter.
	virtual bool import_from_variant( Core::ActionParameterVariant& variant, 
		Core::ActionSource source = Core::ActionSource::NONE_E ) = 0;

	// VALIDATE_VARIANT:
	// Validate a variant parameter
	// This function returns false if the parameter is invalid or cannot be
	// converted and in that case error will describe the error.
	virtual bool validate_variant( Core::ActionParameterVariant& variant, std::string& error ) = 0;

	// -- stateid handling --
public:
	// GET_STATEID:
	// Get the unique id assigned to the state variable
	std::string stateid() const
	{
		return ( stateid_ );
	}

	// GET_BASEID:
	// Get the base id of this state variable
	std::string baseid() const
	{
		return ( stateid_.substr( 0, stateid_.find( ':' ) ) );
	}

	// GET_ID:
	// Get the base id of this state variable
	std::string id() const
	{
		//return ( stateid_.substr( 2, stateid_.find( ":" ) ) );
		return Core::SplitString(  stateid_, "::" )[ 1 ];
	}

protected:
	// INVALIDATE:
	// Invalidate the state variable. This function is called by the "invalidate" function of 
	// StateHandler to release any resource occupied by the state.
	// Default implementation does nothing.
	virtual void invalidate() {}

private:
	std::string stateid_;

	// -- signal handling --
public:
	// STATE_CHANGED_SIGNAL:
	// This signal is triggered when the state is changed
	typedef boost::signals2::signal< void() > state_changed_signal_type;
	state_changed_signal_type state_changed_signal_;

};




} // end namespace Core

#endif

