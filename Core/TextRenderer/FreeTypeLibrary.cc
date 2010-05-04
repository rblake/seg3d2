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

#include <Core/Utils/Exception.h>
#include <Core/Utils/Log.h>
#include <Core/TextRenderer/FreeTypeLibrary.h>

namespace Core
{

FreeTypeLibrary::FreeTypeLibrary( FT_Library library ) :
	library_( library )
{
}

FreeTypeLibrary::~FreeTypeLibrary()
{
	FreeTypeLibraryFactory::DestroyLibrary( this->library_ );
}

FreeTypeFaceHandle FreeTypeLibrary::new_face_from_file( const char* file_name, long index )
{
	FT_Face face;
	FT_Error err = FT_New_Face( this->library_, file_name, index, &face );
	if ( err == 0 )
	{
		return FreeTypeFaceHandle( new FreeTypeFace( face ) );
	}

	SCI_LOG_ERROR( std::string( "Failed to create a new face from font file " ) + file_name );
	return FreeTypeFaceHandle();
}

FreeTypeFaceHandle FreeTypeLibrary::new_face_frome_buffer( const unsigned char* buffer, 
	long size, long index )
{
	FT_Face face;
	FT_Error err = FT_New_Memory_Face( this->library_, buffer, size, index, &face );
	if ( err == 0 )
	{
		return FreeTypeFaceHandle( new FreeTypeFace( face ) );
	}

	SCI_LOG_ERROR( std::string( "Failed to create a new face from buffer " ) );
	return FreeTypeFaceHandle();
}

CORE_SINGLETON_IMPLEMENTATION( FreeTypeLibraryFactory );

FreeTypeLibraryFactory::FreeTypeLibraryFactory()
{
}

FreeTypeLibraryFactory::~FreeTypeLibraryFactory()
{
}

FreeTypeLibraryHandle FreeTypeLibraryFactory::create_library()
{
	lock_type lock( this->get_mutex() );
	FT_Library library;
	FT_Error err = FT_Init_FreeType( &library );
	if ( err == 0 )
	{
		return FreeTypeLibraryHandle( new FreeTypeLibrary( library ) );
	}

	SCI_LOG_ERROR( "Failed to create a FreeType library" );
	return FreeTypeLibraryHandle();
}

void FreeTypeLibraryFactory::destroy_library( FT_Library library )
{
	lock_type lock( this->get_mutex() );
	FT_Done_FreeType( library );
}

FreeTypeLibraryHandle FreeTypeLibraryFactory::CreateLibrary()
{
	return FreeTypeLibraryFactory::Instance()->create_library();
}

void FreeTypeLibraryFactory::DestroyLibrary( FT_Library library )
{
	return  FreeTypeLibraryFactory::Instance()->destroy_library( library );
}

} // end namespace Core