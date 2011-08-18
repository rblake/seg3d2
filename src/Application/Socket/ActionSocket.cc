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

#ifdef _MSC_VER
#pragma warning( disable: 4244 4267 )
#endif

// Boost includes
#include <boost/asio.hpp>
#include <boost/ref.hpp>

// Core includes
#include <Core/Utils/ConnectionHandler.h>
#include <Core/Python/PythonInterpreter.h>

// Application includes
#include <Application/Socket/ActionSocket.h>

namespace Seg3D
{

CORE_SINGLETON_IMPLEMENTATION( ActionSocket );

ActionSocket::ActionSocket() :
	action_socket_thread_( 0 )
{
}

void ActionSocket::start( int portnum )
{
	// Create new thread for running socket code
	action_socket_thread_ = new boost::thread( boost::bind( &ActionSocket::run_action_socket,
	    portnum ) );
}

ActionSocket::~ActionSocket()
{
}

static void WriteOutputToSocket( boost::asio::ip::tcp::socket& socket, std::string output )
{
	boost::asio::write( socket, boost::asio::buffer( output ) );
}

static void WritePromptToSocket( boost::asio::ip::tcp::socket& socket, std::string output )
{
	output = "\r\n" + output;
	WriteOutputToSocket( socket, output );
}

void ActionSocket::run_action_socket( int portnum )
{
	boost::asio::io_service io_service;

	boost::asio::ip::tcp::acceptor acceptor( io_service, boost::asio::ip::tcp::endpoint(
	    boost::asio::ip::tcp::v4(), portnum ) );

	Core::ConnectionHandler connection_handler;

	while ( 1 )
	{
		boost::asio::ip::tcp::socket socket( io_service );
		acceptor.accept( socket );

		// Connect to PythonInterpreter signals
		connection_handler.add_connection( Core::PythonInterpreter::Instance()->prompt_signal_.connect( 
			boost::bind( &WritePromptToSocket, boost::ref( socket ), _1 ) ) );
		connection_handler.add_connection( Core::PythonInterpreter::Instance()->error_signal_.connect( 
			boost::bind( &WriteOutputToSocket, boost::ref( socket ), _1 ) ) );
		connection_handler.add_connection( Core::PythonInterpreter::Instance()->output_signal_.connect( 
			boost::bind( &WriteOutputToSocket, boost::ref( socket ), _1 ) ) );

		boost::asio::write( socket, boost::asio::buffer( std::string( "Welcome to Seg3D\r\n" ) ) );

		boost::system::error_code read_ec;
		while ( !read_ec )
		{
			boost::asio::streambuf buffer;

			boost::asio::read_until( socket, buffer, std::string( "\r\n" ), read_ec );
			std::istream is( &buffer );

			std::string action_string;
			std::getline( is, action_string );

			if ( !read_ec )
			{
				if ( action_string == "exit\r" )
				{
					boost::asio::write( socket, boost::asio::buffer( "Goodbye!\r\n" ) );
					socket.close();
					break;
				}

				Core::PythonInterpreter::Instance()->run_string( action_string );
			}
		}

		connection_handler.disconnect_all();
	}
}

} // end namespace Core

