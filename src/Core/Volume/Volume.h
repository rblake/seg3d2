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

#ifndef CORE_VOLUME_VOLUME_H
#define CORE_VOLUME_VOLUME_H

// Boost includes
#include <boost/utility.hpp>
#include <boost/smart_ptr.hpp>

// Util includes
#include <Core/DataBlock/DataBlock.h>
#include <Core/DataBlock/NrrdData.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Transform.h>
#include <Core/Geometry/GridTransform.h>
#include <Core/Utils/EnumClass.h>

namespace Core
{

// Forward Declaration
class Volume;
typedef boost::shared_ptr< Volume > VolumeHandle;
typedef boost::weak_ptr< Volume > VolumeWeakHandle;

CORE_ENUM_CLASS
(
	VolumeType,
	DATA_E = 0x1, 
	MASK_E = 0x2, 
	LABEL_E = 0x4
)

// Class definition
class Volume : public boost::noncopyable
{
public:
	typedef DataBlock::mutex_type mutex_type;
	typedef DataBlock::lock_type lock_type;
	typedef DataBlock::shared_lock_type shared_lock_type;

	// -- constructor/destructor --
public:

	// Called by the constructors of subclasses
	Volume( const GridTransform& grid_transform );

	virtual ~Volume();

	// -- accessors --
public:

	// GET_TYPE:
	// Get the type of the data layer
	virtual VolumeType get_type() const = 0;

	// GET_GENERATION:
	// Get the generation number of the data block the layer is using
	virtual DataBlock::generation_type get_generation() const = 0;

	// IS_VALID:
	// Check whether the volume has a valid data block
	bool is_valid() const
	{
		return this->get_generation() != -1;
	}

	// GET_GRID_TRANSFORM:
	// Get the grid location
	const GridTransform get_grid_transform() const
	{
		return this->grid_transform_;
	}

	// GET_TRANSFORM:
	// Get the transform of the data
	const Transform get_transform() const
	{
		return this->grid_transform_.transform();
	}
	
	// GET_NX, GET_NY, GET_NZ, GET_SIZE:
	// The volume dimensions
	inline size_t get_nx() const
	{
		return this->nx_;
	}

	inline size_t get_ny() const
	{
		return this->ny_;
	}

	inline size_t get_nz() const
	{
		return this->nz_;
	}

	inline size_t get_size() const
	{
		return this->nx_ * this->ny_ * this->nz_;
	}

	// TO_INDEX:
	// Compute the index into the data blcok based on a coordinate in index space
	inline size_t to_index( size_t x, size_t y, size_t z ) const
	{
		assert( x < this->nx_ && y < this->ny_ && z < this->nz_ );
		return z * this->nxy_ + y * this->nx_ + x;
	}
	
	// APPLY_GRID_TRANSFORM:
	// Compute corresponding point in world coordinates
	Point apply_grid_transform( const Point& pt ) const;

	// APPLY_INVERSE_GRID_TRANSFORM:
	// Compute corresponding point in index coordinates
	Point apply_inverse_grid_transform( const Point& pt ) const;

	// GET_MUTEX:
	// Get the mutex that protects this resource
	virtual mutex_type& get_mutex() = 0;

	// -- internals of volume --
private:

	// Location of the volume object in space
	// NOTE: Currently only axis aligned transforms are allowed
	GridTransform grid_transform_;

	// A precomputed transform in
	Transform inverse_grid_transform_;

	// Cached size information of the GridTransform for fast indexing into the volume
	size_t nx_;
	size_t ny_;
	size_t nz_;
	size_t nxy_;
	
};

} // end namespace Core

#endif
