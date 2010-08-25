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

#ifndef CORE_GRAPHICS_COLORMAP_H
#define CORE_GRAPHICS_COLORMAP_H 

// STL includes
#include <vector>

// Core includes
#include <Core/Geometry/Color.h>
#include <Core/Graphics/Texture.h>

// Boost includes
#include <boost/smart_ptr.hpp> // Needed for shared_ptr

namespace Core 
{

class ColorMap;
typedef boost::shared_ptr< ColorMap > ColorMapHandle;

class ColorMap 
{
public:
	// Creates default rainbow colormap
	ColorMap();

	// Colors should be in range [0, 1]
	ColorMap( const std::vector< Color >& colors, float lookup_min, float lookup_max );

	// TODO Implement copy constructor 

	virtual ~ColorMap();

	// GET_COLORS
	const std::vector< Color >& get_colors() const; 
 
	// SET_COLOR
	// Set/get individual values. Colors should be in range [0, 1].
	void set_color( size_t index, const Color& color );

	// TODO Do we need a function for doing a lookup?  We won't use it for rendering since
	// we'll store the colormap in a texture.  We could include it here anyway for completeness.

	// SET_SIZE
	size_t get_size() const;

	// RESIZE
	// Resize colormap, fill new entries with default color (black)
	void resize( size_t size );

	// SET_LOOKUP_RANGE
	// Set the min/max range of values that will be used to lookup into the ColorMap
	void set_lookup_range( float lookup_min, float lookup_max );

	// GET_LOOKUP_RANGE
	// Get the min/max range of values that will be used to lookup into the ColorMap
	void get_lookup_range( float& lookup_min, float& lookup_max ) const;

	// UPLOAD_TEXTURE
	// Upload colormap data to 1D graphics texture.
	void upload_texture();

	// GET_TEXTURE
	// Get 1D graphics texture with colormap data.
	Texture1DHandle get_texture();

private:
	std::vector< Color > colors_;
	float lookup_min_;
	float lookup_max_;
	Texture1DHandle texture_;
	bool changed_;

	const static Color DEFAULT_COLOR_C;
};

} // End namespace Core


#endif
