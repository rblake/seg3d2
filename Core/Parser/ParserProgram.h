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

#ifndef CORE_PARSER_PARSERPROGRAM_H 
#define CORE_PARSER_PARSERPROGRAM_H 

// STL includes
#include <string>
#include <vector>

// Boost includes
#include <boost/smart_ptr.hpp> // Needed for shared_ptr

// Core includes
#include <Core/Parser/ParserFWD.h>

namespace Core
{

// ParserProgram: This class simply binds a series of expressions together.
class ParserProgram
{
public:
	ParserProgram() {}

	// Add an expression to a program: this is a combination of the raw 
	// unparsed program code and the expression tree
	void add_expression( std::string& expression_string, ParserTreeHandle expression_tree );

	// Retrieve an expression from a program
	// This gets one expression in unparsed and parsed form from the program
	void get_expression( int expression_num, std::string& expression_string,
		ParserTreeHandle& expression_handle );

	void get_expression( int expression_num, ParserTreeHandle& expression_handle );

	// Retrieve the number of expressions in the program
	size_t num_expressions();

	// Add an input variable to the program
	void add_input_variable( std::string name, std::string type = "U", int flags = 0 );

	// Add an output variable to the program
	void add_output_variable( std::string name, std::string type = "U", int flags = 0 );

	void get_input_variables( ParserVariableList& var_list );
	void get_output_variables( ParserVariableList& var_list );

	ParserFunctionCatalogHandle get_catalog();
	void set_catalog( ParserFunctionCatalogHandle catalog );

	// Insert the variables
	void add_const_var( ParserScriptVariableHandle& handle );
	void add_single_var( ParserScriptVariableHandle& handle );
	void add_sequential_var( ParserScriptVariableHandle& handle );

	void add_const_function( ParserScriptFunctionHandle& handle );
	void add_single_function( ParserScriptFunctionHandle& handle );
	void add_sequential_function( ParserScriptFunctionHandle& handle );

	// Get the variables that need to be assigned
	size_t num_const_variables();
	size_t num_single_variables();
	size_t num_sequential_variables();

	bool get_const_variable( size_t j, ParserScriptVariableHandle& handle );
	bool get_single_variable( size_t j, ParserScriptVariableHandle& handle );
	bool get_sequential_variable( size_t j, ParserScriptVariableHandle& handle );

	// Get the functions that need to be assigned
	size_t num_const_functions();
	size_t num_single_functions();
	size_t num_sequential_functions();

	bool get_const_function( size_t j, ParserScriptFunctionHandle& handle );
	bool get_single_function( size_t j, ParserScriptFunctionHandle& handle );
	bool get_sequential_function( size_t j, ParserScriptFunctionHandle& handle );

	// For debugging
	void print();

private:
	// Short cut to the parser function catalog
	ParserFunctionCatalogHandle catalog_;

	// The list of expressions, we store both the raw expression as well as the
	// parsed version in order to effectively report back to user which expression
	// is faulty. As the parsed version is harder to read, the raw version is kept.
	std::vector< std::pair< std::string, ParserTreeHandle > > expressions_;

	// List of variables that exist when the program starts. This contains the
	// list of variables the program can use, without them being defined in one
	// of the expressions
	ParserVariableList input_variables_;

	// List of output variables that need to be generated as they are required
	// by the rest of the program. These variables need to be defined at the
	// end of the program
	ParserVariableList output_variables_;

	// The next series of variables represent the next stage of the parser
	// In this stage everything is a variable or a function, and we have two
	// lists of constants one for double constants and one for string constants

	std::vector< ParserScriptVariableHandle > const_variables_;
	std::vector< ParserScriptVariableHandle > single_variables_;
	std::vector< ParserScriptVariableHandle > sequential_variables_;

	std::vector< ParserScriptFunctionHandle > const_functions_;
	std::vector< ParserScriptFunctionHandle > single_functions_;
	std::vector< ParserScriptFunctionHandle > sequential_functions_;
};

}

#endif
