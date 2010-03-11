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

#ifndef UTILS_GRAPHICS_GLSLPROGRAM_H
#define UTILS_GRAPHICS_GLSLPROGRAM_H

#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#include <GL/glew.h>

#include <Utils/Graphics/GLSLShader.h>

namespace Utils
{

class GLSLProgram;
typedef boost::shared_ptr< GLSLProgram > GLSLProgramHandle;

class GLSLProgram : public boost::noncopyable
{
public:
	GLSLProgram();
	~GLSLProgram();

	void attach_shader( GLSLShaderHandle shader );
	void detach_shader( GLSLShaderHandle shader );

	// Link the program. Returns true if successful, otherwise false.
	// Additional information can be acquired by calling "get_info_log".
	bool link();

	// Validate the program against the current OpenGL state. 
	// Returns true if successful, otherwise false. 
	// Additional information can be acquired by calling "get_info_log".
	bool validate();

	std::string get_info_log();

	void enable();
	void disable();

private:
	GLuint program_id_;
};

} // end namespace Utils

#endif