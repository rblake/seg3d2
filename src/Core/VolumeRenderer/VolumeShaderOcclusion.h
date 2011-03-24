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

#ifndef CORE_VOLUMERENDERER_VOLUMESHADEROCCLUSION_H
#define CORE_VOLUMERENDERER_VOLUMESHADEROCCLUSION_H

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#include <Core/Graphics/ShaderBase.h>

namespace Core
{

class VolumeShaderOcclusion;
typedef boost::shared_ptr< VolumeShaderOcclusion > VolumeShaderOcclusionHandle;

class VolumeShaderOcclusion : public ShaderBase
{
public:
	VolumeShaderOcclusion();
	virtual ~VolumeShaderOcclusion();

	void set_texture_bbox_min( float x, float y, float z );
	void set_texture_bbox_size( float x, float y, float z );
	void set_texel_size( float x, float y, float z );
	void set_voxel_size( float x, float y, float z );
	void set_volume_texture( int tex_unit );
	void set_diffuse_texture( int tex_unit );
	void set_specular_texture( int tex_unit );
	void set_occlusion_sample_texture( int tex_unit );
	void set_occlusion_buffer_texture( int tex_unit );
	void set_num_of_occlusion_samples( int num_of_samples );
	void set_occlusion_extent( float extent );
	void set_slice_distance( float slice_distance );
	void set_clip_plane( const float clip_planes[ 6 ][ 4 ] );
	void set_enable_clip_plane( const int enabled[ 6 ] );
	void set_enable_clipping( bool enabled );

protected:
	virtual bool get_vertex_shader_source( std::string& source );
	virtual bool get_fragment_shader_source( std::string& source );
	virtual bool post_initialize();

private:
	int tex_bbox_min_loc_;
	int tex_bbox_size_loc_;
	int texel_size_loc_;
	int voxel_size_loc_;
	int vol_tex_loc_;
	int diffuse_lut_loc_;
	int specular_lut_loc_;
	int occlusion_sample_lut_loc_;
	int occlusion_buffer_loc_;
	int num_of_occlusion_samples_loc_;
	int occlusion_extent_loc_;
	int slice_distance_loc_;
	int clip_plane_loc_;
	int enable_clip_plane_loc_;
	int enable_clipping_loc_;
};

} // end namespace Seg3D
#endif