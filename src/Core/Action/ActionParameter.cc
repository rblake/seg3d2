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

#include <Core/Action/ActionParameter.h>

namespace Core
{

ActionParameterBase::~ActionParameterBase()
{
}

ActionParameterVariant::ActionParameterVariant()
{
}

ActionParameterVariant::~ActionParameterVariant()
{
}

std::string ActionParameterVariant::export_to_string() const
{
	// Export a value that is still typed or has been convereted to a string
	// if typed_value exist, we need to convert it
	if ( typed_value_.get() )
	{
		return ( typed_value_->export_to_string() );
	}
	else
	{
		// in case typed_value does not exist it must be recorded as a string
		return ( string_value_ );
	}
}

bool ActionParameterVariant::import_from_string( const std::string& str )
{
	// As we do not know the implied type. It can only be recorded as a string
	typed_value_.reset();
	string_value_ = str;

	return ( true );
}

} // namespace Core
