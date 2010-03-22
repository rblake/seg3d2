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

#ifndef UTILS_VOLUME_DATAVOLUMESLICE_H
#define UTILS_VOLUME_DATAVOLUMESLICE_H

#include <Utils/Graphics/PixelBufferObject.h>
#include <Utils/Volume/DataVolume.h>
#include <Utils/Volume/VolumeSlice.h>

namespace Utils
{

class DataVolumeSlice;
typedef boost::shared_ptr< DataVolumeSlice > DataVolumeSliceHandle;

class DataVolumeSlice : public VolumeSlice
{
public:
	DataVolumeSlice( const DataVolumeHandle& data_volume, 
		VolumeSliceType type = VolumeSliceType::AXIAL_E, size_t slice_num = 0 );
	DataVolumeSlice( const DataVolumeSlice& copy );
	virtual ~DataVolumeSlice() {}

	inline double get_data_at ( size_t i, size_t j ) const
	{
		return this->data_block_->get_data_at( this->to_index( i, j ) );
	}

	inline void set_data_at( size_t i, size_t j, double value ) 
	{
		this->data_block_->set_data_at( this->to_index( i, j ), value );
	}

	// Create the texture object
	virtual void initialize_texture();

	// Upload the data slice to graphics texture.
	// NOTE: This function allocates resources on the GPU, so the caller should
	// acquire a lock on the RenderResources before calling this function.
	virtual void upload_texture();

private:
	// Pointer to the data block. The base class keeps a handle of the volume,
	// so it is safe to use a pointer here.
	DataBlock* data_block_;

	// An array of GLenum's for data types, indexed by data_type values
	const static unsigned int GL_DATA_TYPE_C[];

	// An array of GLenum's for GL internal texture formats, indexed by data_type values.
	// These formats are picked to best match the data type.
	const static unsigned int GL_TEXTURE_FORMAT_C[];
};

} // end namespace Utils

#endif