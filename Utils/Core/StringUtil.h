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

#ifndef UTILS_CORE_STRINGUTIL_H
#define UTILS_CORE_STRINGUTIL_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// ConvertString
// This utility contains templated functions for converting strings from
// number and into numbers. The converters handle as well cases missing
// in the standard converter utilities to deal with NaN, inf, and -inf on
// on the input string.


// STL includes
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

namespace Utils
{

// Convert multiple values in a string into a vector with numbers

template< class T >
bool MultipleFromString( const std::string &str, std::vector< T > &values )
{
	values.clear();

	// Clear out any markup of the numbers that make it easier to read and
	// replace it all with spaces.
	std::string data = str;
	for ( size_t j = 0; j < data.size(); j++ )
		if ( ( data[ j ] == '\t' ) || ( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ]
		    == '"' ) || ( data[ j ] == ',' ) || ( data[ j ] == '[' ) || ( data[ j ] == ']' )
		    || ( data[ j ] == '(' ) || ( data[ j ] == ')' ) ) data[ j ] = ' ';

	// Loop over the data and extract all numbers from it.
	for ( size_t p = 0; p < data.size(); )
	{
		// find where the number starts
		while ( ( p < data.size() ) && ( data[ p ] == ' ' ) )
			p++;
		// Exit if we are at the end of the series
		if ( p >= data.size() ) break;

		// strip of the next number
		std::string::size_type next_space = data.find( ' ', p );
		if ( next_space == std::string::npos ) next_space = data.size();

		// Extract the number
		T value;
		if ( FromString( data.substr( p, next_space - p ), value ) ) values.push_back( value );
		p = next_space;

		if ( p >= data.size() ) break;
	}

	// If no numbers were extracted return false
	if ( values.size() > 0 ) return ( true );
	return ( false );
}

// Convert a value into a string

template< class T >
bool FromString( const std::string &str, T &value )
{
	std::string data = str + " ";
	for ( size_t j = 0; j < data.size(); j++ )
		if ( ( data[ j ] == '\t' ) || ( data[ j ] == '\r' ) || ( data[ j ] == '\n' ) || ( data[ j ]
		    == '"' ) || ( data[ j ] == ',' ) || ( data[ j ] == '[' ) || ( data[ j ] == ']' )
		    || ( data[ j ] == '(' ) || ( data[ j ] == ')' ) ) data[ j ] = ' ';

	std::istringstream iss( data );
	iss.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
	try
	{
		iss >> value;
		return ( true );
	}
	catch ( ... )
	{
		return ( false );
	}
}


// Export a value to a string
template< class T >
std::string ToString( T val )
{
	std::ostringstream oss;
	oss << val;
	return ( oss.str() );
}

// Export a value to a string with percision control

template< class T >
std::string ToString( T val, int precision )
{
	std::ostringstream oss;
	oss.precision( precision );
	oss << val;
	return ( oss.str() );
}

// Convert string to upper or lower case 
std::string StringToUpper( std::string );
std::string StringToLower( std::string );

// Special cases that need additional rules to deal with inf and nan
bool FromString( const std::string &str, double &value );
bool FromString( const std::string &str, float &value );

// Functions to strip out spaces at the start or at both ends of the string
void StripSpaces( std::string& str );
void StripSurroundingSpaces( std::string& str );

// Function to split a list of options delimited by a characher into a vector of
// strings
std::vector<std::string> SplitString( const std::string& str, const std::string& delimiter );

} // End namespace Utils

#endif
