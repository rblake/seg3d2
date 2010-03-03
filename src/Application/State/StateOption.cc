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
#include <algorithm>

#include <Application/State/StateOption.h>

namespace Seg3D
{

StateOption::StateOption( const std::string& default_value, const std::string& option_list ) :
	value_( default_value )
{
	// Unwrap the option lists
	std::string option_list_string = Utils::string_to_lower( option_list );
	while ( 1 )
	{
		size_t loc = option_list_string.find( SPLITTER_C );
		if ( loc >= option_list_string.size() )
		{
			option_list_.push_back( option_list_string );
			break;
		}
		option_list_.push_back( option_list_string.substr( 0, loc ) );
		option_list_string = option_list_string.substr( loc + 1 );
	}

	option_list_iterator_type it = 
		std::find( this->option_list_.begin(), this->option_list_.end(), this->value_ );
	if ( it == this->option_list_.end() )
	{
		SCI_THROW_LOGICERROR( std::string( "Option \"" ) + this->value_ +
			"\" not in the option list" );
	}
	this->index_ = static_cast<int>( it - this->option_list_.begin() );
}

StateOption::StateOption( const std::string& default_value,
    const std::vector< std::string >& option_list ) :
	value_( default_value )
{
	option_list_.resize( option_list.size() );
	for ( size_t j = 0; j < option_list.size(); j++ )
	{
		option_list_[ j ] = Utils::string_to_lower( option_list[ j ] );
	}

	option_list_iterator_type it = 
		std::find( this->option_list_.begin(), this->option_list_.end(), this->value_ );
	if ( it == this->option_list_.end() )
	{
		SCI_THROW_LOGICERROR( std::string( "Option \"" ) + this->value_ +
			"\" not in the option list" );
	}
	this->index_ = static_cast<int>( it - this->option_list_.begin() );
}

StateOption::~StateOption()
{
}

std::string StateOption::export_to_string() const
{
	return ( value_ );
}

bool StateOption::set( const std::string& input_value, ActionSource source )
{
	// Lock the state engine so no other thread will be accessing it
	StateEngine::lock_type lock( StateEngine::Instance()->get_mutex() );

	std::string value = Utils::string_to_lower( input_value );
	if ( value != value_ )
	{
		option_list_iterator_type it = 
			std::find( this->option_list_.begin(), this->option_list_.end(), this->value_ );

		if ( option_list_.end() == it )
		{
			if ( source == ACTION_SOURCE_INTERFACE_E )
			{
				// NOTE: This is a special case in which the option requested by the
				// interface does not exist and hence the value may be out of sync and
				// hence needs to be set to the correct value. Hence we generate the
				// signal that indicates that the state has changed anyway.

				// Any other sources are fine as they do not reflect a different value
				// and are validated before the code can reach this point.
				value_changed_signal_( value_, source );
				state_changed_signal_();
			}
			return ( false );
		}
		value_ = value;
		this->index_ = static_cast<int>( it - this->option_list_.begin() );

		value_changed_signal_( value_, source );
		state_changed_signal_();
	}
	return ( true );
}

bool StateOption::import_from_string( const std::string& str, ActionSource source )
{
	std::string value;
	if ( !( Utils::import_from_string( str, value ) ) ) return ( false );

	return ( set( value, source ) );
}

void StateOption::export_to_variant( ActionParameterVariant& variant ) const
{
	variant.set_value( value_ );
}

bool StateOption::import_from_variant( ActionParameterVariant& variant, ActionSource source )
{
	std::string value;
	if ( !( variant.get_value( value ) ) ) return ( false );

	return ( set( value, source ) );
}

bool StateOption::validate_variant( ActionParameterVariant& variant, std::string& error )
{
	std::string value;
	if ( !( variant.get_value( value ) ) )
	{
		error = "Cannot convert the value '" + variant.export_to_string() + "'";
		return ( false );
	}

	value = Utils::string_to_lower( value );
	if ( option_list_.end() == std::find( option_list_.begin(), option_list_.end(), value ) )
	{
		error = "Option '" + value + "' is not a valid option";
		return ( false );
	}

	error = "";
	return ( true );
}

bool StateOption::is_option( const std::string& option )
{
	if ( option_list_.end() == std::find( option_list_.begin(), option_list_.end(),
	    Utils::string_to_lower( option ) ) ) return ( false );
	return ( true );
}

void StateOption::set_option_list( const std::vector< std::string >& option_list )
{
	option_list_ = option_list;
	option_list_iterator_type it = 
		std::find( this->option_list_.begin(), this->option_list_.end(), this->value_ );
	if ( option_list_.end() == it )
	{
		if ( option_list.size() )
		{
			this->index_ = 0;
			value_ = option_list[ 0 ];
		}
		else
		{
			value_ = "";
			this->index_ = -1;
		}
		value_changed_signal_( value_, ACTION_SOURCE_NONE_E );
		state_changed_signal_();
	}
	else
	{
		this->index_ = static_cast<int>( it - this->option_list_.begin() );
	}

	optionlist_changed_signal_();
}

void StateOption::set_option_list( const std::string& option_list )
{
	// Unwrap the option lists
	std::string option_list_string = Utils::string_to_lower( option_list );

	option_list_.clear();
	while ( 1 )
	{
		size_t loc = option_list_string.find( SPLITTER_C );
		if ( loc >= option_list_string.size() )
		{
			option_list_.push_back( option_list_string );
			break;
		}
		option_list_.push_back( option_list_string.substr( 0, loc ) );
		option_list_string = option_list_string.substr( loc + 1 );
	}

	option_list_iterator_type it = 
		std::find( this->option_list_.begin(), this->option_list_.end(), this->value_ );
	if ( option_list_.end() == it )
	{
		if ( option_list.size() ) 
		{
			this->index_ = 0;
			value_ = option_list[ 0 ];
		}
		else 
		{
			value_ = "";
			this->index_ = -1;
		}
		value_changed_signal_( value_, ACTION_SOURCE_NONE_E );
		state_changed_signal_();
	}
	else
	{
		this->index_ = static_cast<int>( it - this->option_list_.begin() );
	}

	optionlist_changed_signal_();
}

void StateOption::set_option_list( const std::string& option_list, const std::string& option )
{
	// Unwrap the option lists
	std::string option_list_string = Utils::string_to_lower( option_list );

	option_list_.clear();
	while ( 1 )
	{
		size_t loc = option_list_string.find( SPLITTER_C );
		if ( loc >= option_list_string.size() )
		{
			option_list_.push_back( option_list_string );
			break;
		}
		option_list_.push_back( option_list_string.substr( 0, loc ) );
		option_list_string = option_list_string.substr( loc + 1 );
	}

	std::string lower_option = Utils::string_to_lower( option );

	option_list_iterator_type it = 
		std::find( this->option_list_.begin(), this->option_list_.end(), lower_option );
	if ( option_list_.end() == it )
	{
		if ( option_list.size() ) 
		{
			this->index_ = 0;
			value_ = option_list[ 0 ];
		}
		else 
		{
			value_ = "";
			this->index_ = -1;
		}
		value_changed_signal_( value_, ACTION_SOURCE_NONE_E );
		state_changed_signal_();
	}
	else
	{
		this->index_ = static_cast<int>( it - this->option_list_.begin() );
		if ( value_ != lower_option )
		{
			value_ = lower_option;
			value_changed_signal_( value_, ACTION_SOURCE_NONE_E );
			state_changed_signal_();
		}
	}

	optionlist_changed_signal_();
}

const char StateOption::SPLITTER_C = '|';

} // end namespace Seg3D
