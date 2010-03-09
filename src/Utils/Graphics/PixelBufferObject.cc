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

#include <Utils/Graphics/PixelBufferObject.h>

namespace Utils
{

PixelBufferObject::PixelBufferObject() :
	BufferObject()
{
}

PixelBufferObject::PixelBufferObject( const BufferObjectHandle& bo ) :
	BufferObject( bo )
{
}

PixelBufferObject::~PixelBufferObject()
{
}

PixelPackBuffer::PixelPackBuffer() :
	PixelBufferObject()
{
	this->target_ = GL_PIXEL_PACK_BUFFER;
	this->query_target_ = GL_PIXEL_PACK_BUFFER_BINDING;
	this->safe_bind();
	this->safe_unbind();
}

PixelPackBuffer::PixelPackBuffer( const BufferObjectHandle& bo ) :
	PixelBufferObject( bo )
{
	this->target_ = GL_PIXEL_PACK_BUFFER;
	this->query_target_ = GL_PIXEL_PACK_BUFFER_BINDING;
}

void PixelPackBuffer::RestoreDefault()
{
	glBindBuffer( GL_PIXEL_PACK_BUFFER, 0 );
}

PixelUnpackBuffer::PixelUnpackBuffer() :
	PixelBufferObject()
{
	this->target_ = GL_PIXEL_UNPACK_BUFFER;
	this->query_target_ = GL_PIXEL_UNPACK_BUFFER_BINDING;
	this->safe_bind();
	this->safe_unbind();
}

PixelUnpackBuffer::PixelUnpackBuffer( const BufferObjectHandle& bo ) :
	PixelBufferObject( bo )
{
	this->target_ = GL_PIXEL_UNPACK_BUFFER;
	this->query_target_ = GL_PIXEL_UNPACK_BUFFER_BINDING;
}

void PixelUnpackBuffer::RestoreDefault()
{
	glBindBuffer( GL_PIXEL_UNPACK_BUFFER, 0 );
}

} // end namespace Utils