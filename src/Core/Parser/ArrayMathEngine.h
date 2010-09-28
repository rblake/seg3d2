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

#ifndef CORE_PARSER_ARRAYMATHENGINE_H 
#define CORE_PARSER_ARRAYMATHENGINE_H 

// Core includes
#include <Core/DataBlock/DataBlock.h>
#include <Core/DataBlock/MaskDataBlock.h>
#include <Core/Parser/ArrayMathInterpreter.h>
#include <Core/Parser/Parser.h>
#include <Core/Parser/ParserFWD.h>

namespace Core
{

class ArrayMathEngine : public Parser, public ArrayMathInterpreter
{
public:
	class OutputDataBlock
	{
	public:
		std::string array_name_;
		std::string data_block_name_; 
		DataBlockHandle data_block_;
	};

	// We don't support outputting a MaskDataBlock because it would have to be registered with the
	// MaskDataBlockManager and locked to prevent conflicts with other masks.  Also, we can avoid 
	// bit operations when copying back the parser result.

public:
	// CALLS TO THIS CLASS SHOULD BE MADE IN THE ORDER
	// THAT THE FUNCTIONS ARE GIVEN HERE

	// Make sure it starts with a clean definition file
	ArrayMathEngine();

	// Generate input arrays
	bool add_input_data_block( std::string name, DataBlockHandle data_block, std::string& error );
	bool add_input_mask_data_block( std::string name, MaskDataBlockHandle mask_data_block, 
		std::string& error );

	//bool add_index( std::string name );
	//bool add_size( std::string name );

	bool add_output_data_block( std::string name, size_t nx, size_t ny, size_t nz, 
		Core::DataType type, std::string& error );
	
	// Setup the expression                        
	bool add_expressions( std::string& expressions );

	bool parse_and_validate( std::string& error );

	// Run the expressions in parallel
	bool run( std::string& error );

	// Extract handles to the results
	bool get_data_block( std::string name, DataBlockHandle& data_block );

	// Clean up the engine
	void clear();

	typedef boost::signals2::signal< void (double) > update_progress_signal_type;

	// UPDATE_PROGRESS:
	// When new information on progress is available this signal is triggered. If this signal is 
	// triggered it should end with a value 1.0 indicating that progress reporting has finised.
	// Progress is measured between 0.0 and 1.0.
	update_progress_signal_type update_progress_signal_;

private:
	void update_progress( double amount );

	// Parser program : the structure of the expressions and simple reduction
	// of the expressions to simple function calls
	ParserProgramHandle pprogram_;
	// Wrapper around the function calls, this piece actually executes the code
	ArrayMathProgramHandle mprogram_;

	// Expression to evaluate before the main expression
	// This one is to extract the variables from the data sources
	// This reduces the amount of actual functions we need to implement
	std::string pre_expression_;
	// The user defined expression
	std::string expression_;
	// Expression to get the data back into the data sinks
	std::string post_expression_;

	// The size of the array engine, the first call that add an array that is
	// bigger than 1, will set this variable
	// Any subsequent array that does not match the size will cause an error
	size_type array_size_;

	// Data that needs to be stored as it is needed before and after the parser is
	// done. An output needs to be set before the parser, otherwise it is optimized
	// away, but the type is only know when the parser has validated and optimized
	// the expression tree
	std::vector< OutputDataBlock > data_block_data_;
};

}

#endif
