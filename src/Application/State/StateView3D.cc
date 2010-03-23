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

// Application includes
#include <Application/State/StateView3D.h>

namespace Seg3D
{

StateView3D::StateView3D( const std::string& stateid ) :
	StateViewBase( stateid )
{
}

StateView3D::~StateView3D()
{
}

std::string StateView3D::export_to_string() const
{
	return ( Utils::ExportToString( value_ ) );
}

bool StateView3D::import_from_string( const std::string& str, ActionSource source )
{
	Utils::View3D value;
	if ( !( Utils::ImportFromString( str, value ) ) ) return ( false );

	return ( set( value, source ) );
}

bool StateView3D::set( const Utils::View3D& value, ActionSource source )
{
	// Lock the state engine so no other thread will be accessing it
	StateEngine::lock_type lock( StateEngine::Instance()->get_mutex() );

	if ( value != value_ )
	{
		value_ = value;
		value_changed_signal_( value_, source );
		state_changed_signal_();
	}
	return ( true );
}

void StateView3D::rotate( const Utils::Vector& axis, double angle )
{
	{
		// Lock the state engine so no other thread will be accessing it
		StateEngine::lock_type lock( StateEngine::Instance()->get_mutex() );
		this->value_.rotate( axis, angle );
	}
	this->state_changed_signal_();
}

void StateView3D::scale( double ratio )
{
	{
		// Lock the state engine so no other thread will be accessing it
		StateEngine::lock_type lock( StateEngine::Instance()->get_mutex() );
		this->value_.scale( ratio );
	}
	this->state_changed_signal_();
}

void StateView3D::translate( const Utils::Vector& offset )
{
	{
		// Lock the state engine so no other thread will be accessing it
		StateEngine::lock_type lock( StateEngine::Instance()->get_mutex() );
		this->value_.translate( offset );
	}
	this->state_changed_signal_();
}

void StateView3D::export_to_variant( ActionParameterVariant& variant ) const
{
	variant.set_value( value_ );
}

bool StateView3D::import_from_variant( ActionParameterVariant& variant, ActionSource source )
{
	Utils::View3D value;
	if ( !( variant.get_value( value ) ) ) return ( false );

	return ( set( value, source ) );
}

bool StateView3D::validate_variant( ActionParameterVariant& variant, std::string& error )
{
	Utils::View3D value;
	if ( !( variant.get_value( value ) ) )
	{
		error = "Cannot convert the value '" + variant.export_to_string()
		    + "' to a 3D Camera position";
		return ( false );
	}

	error = "";
	return ( true );
}

} // end namespace Seg3D
