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
#pragma warning( disable: 4244 )
#endif

#include <Python.h>
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/thread/mutex.hpp>

#include <Core/Utils/Exception.h>
#include <Core/Utils/Lockable.h>

#include <Application/PythonModule/ActionModule.h>
#include <Application/PythonModule/PythonInterpreter.h>

//////////////////////////////////////////////////////////////////////////
// Python flags defined in pythonrun.c
//////////////////////////////////////////////////////////////////////////

// Py_IgnoreEnvironmentFlag:
// Set this flag to tell python to ignore environment variables.
extern int Py_IgnoreEnvironmentFlag;

// Py_InteractiveFlag:
// Set this flag to run the interpreter in interactive mode when Py_Main is called.
extern int Py_InteractiveFlag;

// Py_InspectFlag:
// Set this flag so the program won't exit at SystemError
extern int Py_InspectFlag;

// Py_NoSiteFlag:
// Suppress 'import site'
extern int Py_NoSiteFlag;

// Py_OptimizeFlag:
extern int Py_OptimizeFlag;

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// Class PythonInterpreterPrivate
//////////////////////////////////////////////////////////////////////////

class PythonInterpreterPrivate : public Core::Lockable
{
public:
	std::string read_from_console( const int bytes = -1 );

	// The name of the executable
	wchar_t* program_name_;
	// An instance of python CommandCompiler object (defined in codeop.py)
	boost::python::object compiler_;
	// The context of the Python main module.
	boost::python::object globals_;
	// Whether the Python interpreter has been initialized.
	bool initialized_;
	bool terminal_running_;
	// The action context for the Python thread.
	PythonActionContextHandle action_context_;
	// The input buffer
	std::string input_buffer_;
	// Whether the interpreter is waiting for input
	bool waiting_for_input_;
	// The command buffer (for Python statement split into multiple lines)
	std::string command_buffer_;
	// Python sys.ps1
	std::string prompt1_;
	// Python sys.ps2
	std::string prompt2_;

	// Condition variable to make sure the PythonInterpreter thread has 
	// completed initialization before continuing the main thread.
	boost::condition_variable thread_condition_variable_;
};

std::string PythonInterpreterPrivate::read_from_console( const int bytes /*= -1 */ )
{
	lock_type lock( this->get_mutex() );

	std::string result;
	if ( !input_buffer_.empty() )
	{
		if ( bytes <= 0 )
		{
			result = input_buffer_;
			input_buffer_.clear();
		}
		else
		{
			result = input_buffer_.substr( 0, bytes );
			if ( bytes < input_buffer_.size() )
			{
				input_buffer_ = input_buffer_.substr( bytes );
			}
			else
			{
				input_buffer_.clear();
			}
		}
	}

	int more_bytes = bytes - static_cast< int >( result.size() );
	while ( ( bytes <= 0 && result.empty() ) ||
		( bytes > 0 && more_bytes > 0 ) )
	{
		this->waiting_for_input_ = true;
		this->thread_condition_variable_.wait( lock );

		// Abort reading if an interrupt signal has been received.
		if ( PyErr_CheckSignals() != 0 ) break;

		if ( bytes <= 0 )
		{
			result = input_buffer_;
			input_buffer_.clear();
		}
		else
		{
			result += input_buffer_.substr( 0, more_bytes );
			if ( more_bytes < input_buffer_.size() )
			{
				input_buffer_ = input_buffer_.substr( more_bytes );
			}
			else
			{
				input_buffer_.clear();
			}
		}

		more_bytes = bytes - static_cast< int >( result.size() );
	}

	this->waiting_for_input_ = false;
	return result;
}

} // end namespace Seg3D


//////////////////////////////////////////////////////////////////////////
// Class PythonTerminal
//////////////////////////////////////////////////////////////////////////
class PythonStdIO
{
public:
	boost::python::object read( int n )
	{
		std::string data = Seg3D::PythonInterpreter::Instance()->private_->read_from_console( n );
		boost::python::str pystr( data.c_str() );
		return pystr.encode();
	}

	std::string readline()
	{
		return Seg3D::PythonInterpreter::Instance()->private_->read_from_console();
	}

	int write( std::string data )
	{
		Seg3D::PythonInterpreter::Instance()->output_signal_( data );
		return static_cast< int >( data.size() );
	}
};

class PythonStdErr
{
public:
	int write( std::string data )
	{
		Seg3D::PythonInterpreter::Instance()->error_signal_( data );
		return static_cast< int >( data.size() );
	}
};

BOOST_PYTHON_MODULE( interpreter )
{
	boost::python::class_< PythonStdIO >( "terminalio" )
		.def( "read", &PythonStdIO::read )
		.def( "readline", &PythonStdIO::readline )
		.def( "write", &PythonStdIO::write );

	boost::python::class_< PythonStdErr >( "terminalerr" )
		.def( "write", &PythonStdErr::write );
}

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// Class PythonInterpreter
//////////////////////////////////////////////////////////////////////////

CORE_SINGLETON_IMPLEMENTATION( PythonInterpreter );
	
PythonInterpreter::PythonInterpreter() :
	Core::EventHandler(),
	private_( new PythonInterpreterPrivate )
{
	this->private_->program_name_ = L"";
	this->private_->initialized_ = false;
	this->private_->terminal_running_ = false;
	this->private_->waiting_for_input_ = false;
	this->private_->action_context_.reset( new PythonActionContext );
}

PythonInterpreter::~PythonInterpreter()
{
	// NOTE: Boost.Python requires that we don't call Py_Finalize
	//Py_Finalize();
}

void PythonInterpreter::initialize_eventhandler()
{
	using namespace boost::python;

	PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
	PyImport_AppendInittab( "seg3d", PyInit_seg3d );
	PyImport_AppendInittab( "interpreter", PyInit_interpreter );
	Py_SetProgramName( this->private_->program_name_ );
	boost::filesystem::path lib_path( this->private_->program_name_ );
	lib_path = lib_path.parent_path() / PYTHONPATH;
	Py_SetPath( lib_path.wstring().c_str() );
	Py_IgnoreEnvironmentFlag = 1;
	Py_InspectFlag = 1;
	Py_OptimizeFlag = 2;
	Py_NoSiteFlag = 1;
	//Py_InteractiveFlag = 1;
	Py_Initialize();

	// Import seg3d module and everything inside
	PyRun_SimpleString( "import seg3d\n"
		"from seg3d import *\n" );

	// Create the compiler object
	PyRun_SimpleString( "from codeop import CommandCompiler\n"
		"__internal_compiler = CommandCompiler()\n" );
 	boost::python::object main_module = boost::python::import( "__main__" );
 	boost::python::object main_namespace = main_module.attr( "__dict__" );
	this->private_->compiler_ = main_namespace[ "__internal_compiler" ];
	this->private_->globals_ = main_namespace;

	// Set up the prompt strings
	PyRun_SimpleString( "import sys\n"
		"try:\n"
		"\tsys.ps1\n"
		"except AttributeError:\n"
		"\tsys.ps1 = \">>> \"\n"
		"try:\n"
		"\tsys.ps2\n"
		"except AttributeError:\n"
		"\tsys.ps2 = \"... \"\n" );
	boost::python::object sys_module = main_namespace[ "sys" ];
	boost::python::object sys_namespace = sys_module.attr( "__dict__" );
	this->private_->prompt1_ = boost::python::extract< std::string >( sys_namespace[ "ps1" ] );
	this->private_->prompt2_ = boost::python::extract< std::string >( sys_namespace[ "ps2" ] );

	// Hook up the I/O
	PyRun_SimpleString( "import interpreter\n"
		"__term_io = interpreter.terminalio()\n"
		"__term_err = interpreter.terminalerr()\n" );
	PyRun_SimpleString( "import sys\n" 
		"sys.stdin = __term_io\n"
		"sys.stdout = __term_io\n"
		"sys.stderr = __term_err\n" );

	// Remove intermediate python variables
	PyRun_SimpleString( "del (interpreter, __internal_compiler, __term_io, __term_err)\n" );

	this->private_->thread_condition_variable_.notify_one();
}

void PythonInterpreter::initialize( wchar_t* program_name )
{
	this->private_->program_name_ = program_name;

	PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
	this->start_eventhandler();
	this->private_->thread_condition_variable_.wait( lock );
	this->private_->initialized_ = true;
}

void PythonInterpreter::print_banner()
{
	if ( !this->is_eventhandler_thread() )
	{
		this->post_event( boost::bind( &PythonInterpreter::print_banner, this ) );
		return;
	}

	PyRun_SimpleString( "print('Pyton %s on %s' % (sys.version, sys.platform))\n" );
	this->prompt_signal_( this->private_->prompt1_ );
}

void PythonInterpreter::run_string( std::string command )
{
	{
		PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
		if ( !this->private_->initialized_ )
		{
			CORE_THROW_LOGICERROR( "The python interpreter hasn't been initialized!" );
		}
	}
	
	if ( !this->is_eventhandler_thread() )
	{
		{
			PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
			// If the Python thread is currently waiting for input, feed the string
			// to the input buffer directly and return.
			if ( this->private_->waiting_for_input_ )
			{
				this->private_->input_buffer_ = command + "\n";
				this->private_->thread_condition_variable_.notify_one();
				return;
			}
		}

		this->post_event( boost::bind( &PythonInterpreter::run_string, this, command ) );
		return;
	}

	// Clear any previous Python errors.
	PyErr_Clear();

	// Append the command to the current command buffer
	if ( this->private_->command_buffer_.empty() )
	{
		this->private_->command_buffer_ = command;
	}
	else
	{
		if ( *this->private_->command_buffer_.rbegin() != '\n' )
		{
			this->private_->command_buffer_ += "\n";
		}
		this->private_->command_buffer_ += command;
	}

	// Compile the statement in the buffer
	boost::python::object code_obj;
	try
	{
		code_obj = this->private_->compiler_( this->private_->command_buffer_ );
	}
	catch ( ... ) {}

	// If an error happened during compilation, print the error message
	if ( PyErr_Occurred() != NULL )
	{
		PyErr_Print();
	}
	// If compilation succeeded and the code object is not Py_None
	else if ( code_obj )
	{
		PyObject* result = PyEval_EvalCode( code_obj.ptr(), this->private_->globals_.ptr(), NULL );
		Py_XDECREF( result );
		if ( PyErr_Occurred() != NULL )
		{
			if ( PyErr_ExceptionMatches( PyExc_EOFError ) )
			{
				this->error_signal_( "\nKeyboardInterrupt\n" );
				PyErr_Clear();
			}
			else
			{
				PyErr_Print();
			}
		}
	}
	// If the code object is Py_None, prompt for more input
	else
	{
		this->prompt_signal_( this->private_->prompt2_ );
		return;
	}

	this->private_->command_buffer_.clear();
	this->prompt_signal_( this->private_->prompt1_ );
}

void PythonInterpreter::run_file( std::string file_name )
{
	{
		PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
		if ( !this->private_->initialized_ )
		{
			CORE_THROW_LOGICERROR( "The python interpreter hasn't been initialized!" );
		}
	}

	if ( !this->is_eventhandler_thread() )
	{
		this->post_event( boost::bind( &PythonInterpreter::run_file, this, file_name ) );
		return;
	}

	FILE* fp = fopen( file_name.c_str(), "r" );
	if ( fp != 0 )
	{
		PyRun_SimpleFileEx( fp, file_name.c_str(), 1 );
	}
}

void PythonInterpreter::interrupt()
{
	if ( !this->is_eventhandler_thread() )
	{
		PyErr_SetInterrupt();
		this->private_->thread_condition_variable_.notify_all();
		this->post_event( boost::bind( &PythonInterpreter::interrupt, this ) );
	}
	else
	{
		if ( PyErr_CheckSignals() != 0 )
		{
			this->error_signal_( "\nKeyboardInterrupt\n" );
			this->private_->command_buffer_.clear();
			this->prompt_signal_( this->private_->prompt1_ );
		}
	}
}

void PythonInterpreter::start_terminal()
{
	{
		PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
		if ( !this->private_->initialized_ )
		{
			CORE_THROW_LOGICERROR( "The python interpreter hasn't been initialized!" );
		}

		if ( this->private_->terminal_running_ )
		{
			return;
		}
	}

	if ( !this->is_eventhandler_thread() )
	{
		this->post_event( boost::bind( &PythonInterpreter::start_terminal, this ) );
		return;
	}
	
	PythonInterpreterPrivate::lock_type lock( this->private_->get_mutex() );
	this->private_->terminal_running_ = true;
	lock.unlock();

	//PyRun_SimpleString("import sys\n");
	//PyRun_SimpleString("import io\n");
	//PyRun_SimpleString("sys.stdout=io.open('CONOUT$', 'wt')\n");
	//PyObject* sys = PyImport_ImportModule("sys");
	//PyObject* io = PyImport_ImportModule("io");
	//PyObject* pystdout = PyObject_CallMethod(io, "open", "ss", "CONOUT$", "wt");
	//if (-1 == PyObject_SetAttrString(sys, "stdout", pystdout)) {
	//	/* Announce your error to the world */
	//}
	//Py_DECREF(sys);
	//Py_DECREF(io);
	//Py_DECREF(pystdout);
	wchar_t** argv = new wchar_t*[ 2 ];
	argv[ 0 ] = this->private_->program_name_;
	argv[ 1 ] = 0;
	Py_Main( 1, argv );
	delete[] argv;
	
	lock.lock();
	this->private_->terminal_running_ = false;
}

PythonActionContextHandle PythonInterpreter::get_action_context()
{
	return this->private_->action_context_;
}

PythonActionContextHandle PythonInterpreter::GetActionContext()
{
	return Instance()->get_action_context();
}

} // end namespace Core
