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

#include <Utils/Core/Log.h>

#include <Application/Action/ActionContext.h>

namespace Seg3D
{

ActionContext::ActionContext()
{
}

ActionContext::~ActionContext()
{
}

void ActionContext::report_error( const std::string& error )
{
	SCI_LOG_ERROR(error);
}

void ActionContext::report_warning( const std::string& warning )
{
	SCI_LOG_WARNING(warning);
}

void ActionContext::report_message( const std::string& message )
{
	SCI_LOG_MESSAGE(message);
}

void ActionContext::report_result( const ActionResultHandle& result )
{
}

void ActionContext::report_status( ActionStatus status )
{
	status_ = status;
}

void ActionContext::report_need_resource( ResourceLockHandle& resource )
{
}

ActionStatus ActionContext::status()
{
	return ( status_ );
}

ActionSource ActionContext::source()
{
	return ( ACTION_SOURCE_COMMANDLINE_E );
}

void ActionContext::report_done()
{
}

} // end namespace Seg3D
