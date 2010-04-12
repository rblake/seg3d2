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

#include <cassert>

#include <Utils/Core/Exception.h>
#include <Utils/Core/Log.h>
#include <Utils/TextRenderer/FreeTypeGlyph.h>

namespace Utils
{

FreeTypeGlyph::FreeTypeGlyph( FT_Glyph glyph ) :
	glyph_( glyph )
{
}

FreeTypeGlyph::FreeTypeGlyph( const FreeTypeGlyph& copy )
{
	FT_Glyph_Copy( copy.glyph_, &this->glyph_ );
}

FreeTypeGlyph::~FreeTypeGlyph()
{
	FT_Done_Glyph( this->glyph_ );
}

void FreeTypeGlyph::transform( FT_Matrix* matrix, FT_Vector* delta )
{
	FT_Glyph_Transform( this->glyph_, matrix, delta );
}

FreeTypeBitmapGlyphHandle FreeTypeGlyph::render_to_bitmap( FT_Vector* origin )
{
	FT_Glyph glyph;

	if ( this->glyph_->format == FT_GLYPH_FORMAT_BITMAP )
	{
		FT_Glyph_Copy( this->glyph_, &glyph );
	}
	else
	{
		glyph = this->glyph_;
		FT_Glyph_To_Bitmap( &glyph, FT_RENDER_MODE_NORMAL, origin, 0 );
	}

	return FreeTypeBitmapGlyphHandle( new FreeTypeBitmapGlyph( ( FT_BitmapGlyph )glyph ) );
}

FreeTypeBitmapGlyph::FreeTypeBitmapGlyph( FT_BitmapGlyph glyph_bitmap ) :
	glyph_bitmap_( glyph_bitmap )
{
}

FreeTypeBitmapGlyph::~FreeTypeBitmapGlyph()
{
	FT_Done_Glyph( (FT_Glyph)this->glyph_bitmap_ );
}

void FreeTypeBitmapGlyph::draw( unsigned char* target, const int width, const int height, 
	const int x_offset, const int y_offset, const float red, const float green, 
	const float blue, const float alpha, bool blend )
{
	FT_Bitmap& bitmap = this->glyph_bitmap_->bitmap;
	assert( bitmap.pixel_mode == FT_PIXEL_MODE_GRAY );

	int pen_y = y_offset + this->glyph_bitmap_->top;
	unsigned char* pixels = bitmap.buffer;
	for ( int i = 0; i < bitmap.rows; i++, pen_y-- )
	{
		if ( pen_y >= height ) continue;
		if ( pen_y < 0 ) break;

		int pen_x = x_offset + this->glyph_bitmap_->left;
		for ( int j = 0; j < bitmap.width; j++, pen_x++ )
		{
			if ( pen_x < 0 ) continue;
			if ( pen_x >= width ) break;

			float gray_scale = pixels[ j ] / float( bitmap.num_grays - 1 );
			float src_alpha = alpha * gray_scale;
			int pixel_pos = ( pen_y * width + pen_x ) * 4;
			if ( blend )
			{
				float dst_factor = 1 - src_alpha;
//				target[ pixel_pos ] = static_cast< unsigned char >( 
//					( red * src_alpha + target[ pixel_pos ] / 255.0 * dst_factor ) * 255 );
//				target[ pixel_pos + 1 ] = static_cast< unsigned char >( 
//					( green * src_alpha + target[ pixel_pos + 1 ] / 255.0 * dst_factor ) * 255 );
//				target[ pixel_pos + 2 ] = static_cast< unsigned char >( 
//					( blue * src_alpha + target[ pixel_pos + 2 ] / 255.0 * dst_factor ) * 255 );
//				target[ pixel_pos + 3 ] = static_cast< unsigned char >( 
//					( 1 - dst_factor * ( 1 - target[ pixel_pos + 3 ] / 255.0 ) ) * 255 );
				target[ pixel_pos ] = static_cast< unsigned char >( red * 255 );
				target[ pixel_pos + 1 ] = static_cast< unsigned char >( green * 255 );
				target[ pixel_pos + 2 ] = static_cast< unsigned char >( blue * 255 );
				target[ pixel_pos + 3 ] = static_cast< unsigned char >( 
					( gray_scale + dst_factor * target[ pixel_pos + 3 ] / 255.0 ) * 255 );
			}
			else
			{
				target[ pixel_pos ] = static_cast< unsigned char >( red * 255 );
				target[ pixel_pos + 1 ] = static_cast< unsigned char >( green * 255 );
				target[ pixel_pos + 2 ] = static_cast< unsigned char >( blue * 255 );
				target[ pixel_pos + 3 ] = static_cast< unsigned char >( src_alpha * 255 );
			}
		}
		pixels += bitmap.pitch;
	}
}

} // end namespace Utils