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

#ifndef CORE_ISOSURFACE_ISOSURFACEBASE_H
#define CORE_ISOSURFACE_ISOSURFACEBASE_H

// STL includes
#include <vector>

// Boost includes
#include <boost/smart_ptr.hpp> // Needed for shared_ptr
#include <boost/utility.hpp> // Needed for noncopyable

// Core includes
#include <Core/Graphics/ColorMap.h>
#include <Core/Utils/Lockable.h>

namespace Core
{

class IsosurfaceBase;
typedef boost::shared_ptr< IsosurfaceBase > IsosurfaceBaseHandle;

class IsosurfaceBase : public Core::RecursiveLockable
{
public:
  // COMPUTE:
  /// Compute isosurface.  quality_factor must be one of: {0.125, 0.25, 0.5, 1.0}
  virtual void compute( double quality_factor, bool capping_enabled, boost::function< bool () > check_abort ) = 0;

  // SURFACE_AREA:
  /// Return the area of the isosurface.
  virtual double surface_area() const = 0;

  // SET_COLOR_MAP:
  /// Set mapping from vertex values to RGB colors.
  /// NOTE: This function is not thread-safe. Passing handle since colormap is unlikely
  /// to be modified after creation.
  virtual void set_color_map( ColorMapHandle color_map ) = 0;

  // GET_COLOR_MAP:
  /// Get mapping from vertex values to RGB colors
  virtual ColorMapHandle get_color_map() const = 0;

  // REDRAW:
  /// Render the isosurface.  This function doesn't work in isolation -- it must be called from the
  /// Seg3D Renderer.
  virtual void redraw( bool use_colormap ) = 0;

	// EXPORT_LEGACY_ISOSURFACE:
	/// Write points to .pts file, faces to .fac file, and values (if assigned) to .val file.
	/// Returns true on success, false on failure.
	///
	/// path: Path to existing directory where files should be written.
	/// file_prefix: File prefix to use for output files (no extension).
	///
	/// Format for .pts:
	/// x y z
	/// x y z
	/// ...
	///
	/// Format for .fac:
	/// p1 p2 p3
	/// p1 p2 p3
	/// ...
	///
	/// Format for .val:
	/// v1
	/// v2
	/// ...
	///
	/// Note: can't call this function "export" because it is reserved by the Visual C++ compiler.
	virtual bool export_legacy_isosurface( const boost::filesystem::path& path, const std::string& file_prefix ) = 0;

  // EXPORT_VTK_ISOSURFACE:
  /// Writes out an isosurface in VTK mesh format
  virtual bool export_vtk_isosurface( const boost::filesystem::path& filename ) = 0;

  // EXPORT_STL_ISOSURFACE:
  /// Writes out an isosurface in STL file format
  virtual bool export_stl_isosurface( const boost::filesystem::path& filename, const std::string& name ) = 0;
};

} // end namespace Core

#endif