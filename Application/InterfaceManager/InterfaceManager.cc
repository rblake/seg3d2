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

#include <Core/Application/Application.h>
#include <Core/Interface/Interface.h>

#include <Application/InterfaceManager/InterfaceManager.h>

namespace Seg3D
{

const size_t InterfaceManager::version_number_ = 1;

CORE_SINGLETON_IMPLEMENTATION( InterfaceManager );

InterfaceManager::InterfaceManager() :
	StateHandler( "interface", version_number_, false )
{
	// set up state variables
	add_state( "fullscreen", full_screen_state_, false );
}

InterfaceManager::~InterfaceManager()
{
	disconnect_all();
}

void InterfaceManager::add_windowid( const std::string& windowid )
{
	std::string lower_windowid = Core::StringToLower( windowid );
	boost::unique_lock< boost::mutex > lock( windowid_list_mutex_ );
	if ( windowid_list_.find( lower_windowid ) == windowid_list_.end() )
	{
		windowid_list_.insert( lower_windowid );
	}
}

bool InterfaceManager::is_windowid( const std::string& windowid )
{
	std::string lower_windowid = Core::StringToLower( windowid );
	boost::unique_lock< boost::mutex > lock( windowid_list_mutex_ );
	return ( windowid_list_.find( lower_windowid ) != windowid_list_.end() );
}

} // end namespace Seg3D
