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

#ifndef CORE_ACTION_ACTIONPARAMETERTYPE_H
#define CORE_ACTION_ACTIONPARAMETERTYPE_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// STL
#include <string>

// Core includes
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

namespace Core
{

// PARAMETERTYPE:
// Function that returns the type of the parameter

const int PARAMETER_UNKNOWN_C = -1;
const int PARAMETER_CHAR_C = 1;
const int PARAMETER_UCHAR_C = 2;
const int PARAMETER_SHORT_C = 3;
const int PARAMETER_USHORT_C = 4;
const int PARAMETER_INT_C = 5;
const int PARAMETER_UINT_C = 6;
const int PARAMETER_LONGLONG_C = 7;
const int PARAMETER_ULONGLONG_C = 8;
const int PARAMETER_FLOAT_C = 9;
const int PARAMETER_DOUBLE_C = 10;
const int PARAMETER_STRING_C = 11;
const int PARAMETER_POINT_C = 12;
const int PARAMETER_VECTOR_C = 13;

template< class T >
inline int ParameterType( T& type ) {  return PARAMETER_UNKNOWN_C; }

inline int ParameterType( char& type ) { return PARAMETER_CHAR_C; }
inline int ParameterType( unsigned char& type ) { return PARAMETER_UCHAR_C; }
inline int ParameterType( short& type ) { return PARAMETER_SHORT_C; }
inline int ParameterType( unsigned short& type ) { return PARAMETER_USHORT_C; }
inline int ParameterType( int& type ) { return PARAMETER_INT_C; }
inline int ParameterType( unsigned int& type ) { return PARAMETER_UINT_C; }
inline int ParameterType( long long& type ) { return PARAMETER_LONGLONG_C; }
inline int ParameterType( unsigned long long& type ) { return PARAMETER_ULONGLONG_C; }
inline int ParameterType( float& type ) { return PARAMETER_FLOAT_C; }
inline int ParameterType( double& type ) { return PARAMETER_DOUBLE_C; }
inline int ParameterType( std::string& type ) { return PARAMETER_STRING_C; }
inline int ParameterType( Point& type ) { return PARAMETER_POINT_C; }
inline int ParameterType( Vector& type ) { return PARAMETER_VECTOR_C; }

} // namespace Core

#endif
