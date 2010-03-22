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

#include <Utils/Math/MathFunctions.h>
#include <Utils/Volume/VolumeSlice.h>

namespace Utils
{

VolumeSlice::VolumeSlice( const VolumeHandle& volume, 
						 VolumeSliceType type, size_t slice_num ) :
	slice_changed_( true ), 
	size_changed_( true ),
	volume_( volume ), 
	slice_type_( type ), 
	slice_number_ ( slice_num )
{
	this->update_position();
	this->slice_number_ = Min( this->slice_number_, this->number_of_slices_ - 1 );
}

VolumeSlice::VolumeSlice( const VolumeSlice& copy ) :
	slice_changed_( copy.slice_changed_ ),
	size_changed_( copy.size_changed_ ),
	nx_( copy.nx_ ),
	ny_( copy.ny_ ),
	number_of_slices_( copy.number_of_slices_ ),
	left_( copy.left_ ), right_( copy.right_ ), 
	bottom_( copy.bottom_ ), top_( copy.top_ ),
	bottom_left_( copy.bottom_left_ ),
	bottom_right_( copy.bottom_right_ ),
	top_left_( copy.top_left_ ),
	top_right_( copy.top_right_ ),
	texture_( copy.texture_ ),
	volume_( copy.volume_ ),
	slice_type_( copy.slice_type_ ),
	slice_number_( copy.slice_number_ )
{
}

VolumeSlice::~VolumeSlice()
{
	this->disconnect_all();
}

void VolumeSlice::set_slice_type( VolumeSliceType type )
{
	if ( this->slice_type_ != type )
	{
		this->slice_changed_ = true;
		this->size_changed_ = true;
		this->slice_type_ = type;

		this->update_position();
		this->slice_number_ = Min( this->slice_number_, this->number_of_slices_ - 1 );
	}
}

void VolumeSlice::set_slice_number( size_t slice_num )
{
	slice_num = Min( slice_num, this->number_of_slices_ - 1 );
	if ( this->slice_number_ != slice_num )
	{
		this->slice_number_ = slice_num;
		this->slice_changed_ = true;
	}
}

void VolumeSlice::update_position()
{
	Point index;
	this->to_index( 0, 0, index );
	this->bottom_left_ = this->volume_->apply_grid_transform( index );
	this->to_index( this->nx_ - 1, 0, index );
	this->bottom_right_ = this->volume_->apply_grid_transform( index );
	this->to_index( this->nx_ - 1, this->ny_ - 1, index );
	this->top_right_ = this->volume_->apply_grid_transform( index );
	this->to_index( 0, this->ny_ - 1, index );
	this->top_left_ = this->volume_->apply_grid_transform( index );

	switch( this->slice_type_ )
	{
	case VolumeSliceType::AXIAL_E:
		this->nx_ = this->volume_->nx();
		this->ny_ = this->volume_->ny();
		this->number_of_slices_ = this->volume_->nz();
		this->left_ = this->bottom_left_.x();
		this->right_ = this->top_right_.x();
		this->bottom_ = this->bottom_left_.y();
		this->top_ = this->top_right_.y();
		break;
	case VolumeSliceType::CORONAL_E:
		this->nx_ = this->volume_->nz();
		this->ny_ = this->volume_->nx();
		this->number_of_slices_ = this->volume_->ny();
		this->left_ = this->bottom_left_.z();
		this->right_ = this->top_right_.z();
		this->bottom_ = this->bottom_left_.x();
		this->top_ = this->top_right_.x();
		break;
	case VolumeSliceType::SAGITTAL_E:
		this->nx_ = this->volume_->ny();
		this->ny_ = this->volume_->nz();
		this->number_of_slices_ = this->volume_->nx();
		this->left_ = this->bottom_left_.y();
		this->right_ = this->top_right_.y();
		this->bottom_ = this->bottom_left_.z();
		this->top_ = this->top_right_.z();
		break;
	default:
		assert( false );
		break;
	}

	// TODO: remove this. It's for testing only
	this->slice_number_ = this->number_of_slices_ / 2;
}

size_t VolumeSlice::to_index( size_t i, size_t j ) const
{
	switch ( this->slice_type_ )
	{
	case VolumeSliceType::AXIAL_E:
		return this->volume_->to_index( i, j, this->slice_number_ );
	case VolumeSliceType::CORONAL_E:
		return this->volume_->to_index( j, this->slice_number_, i );
	default:
		return this->volume_->to_index( this->slice_number_, i, j );
	}
}

void VolumeSlice::to_index( size_t i, size_t j, Point& index ) const
{
	switch ( this->slice_type_ )
	{
	case VolumeSliceType::AXIAL_E:
		index[ 0 ] = static_cast<double>( i );
		index[ 1 ] = static_cast<double>( j );
		index[ 2 ] = static_cast<double>( this->slice_number_ );
		break;
	case VolumeSliceType::CORONAL_E:
		index[ 0 ] = static_cast<double>( j );
		index[ 1 ] = static_cast<double>( this->slice_number_ );
		index[ 2 ] = static_cast<double>( i );
		break;
	default:
		index[ 0 ] = static_cast<double>( this->slice_number_ );
		index[ 1 ] = static_cast<double>( i );
		index[ 2 ] = static_cast<double>( j );
	}
}

} // end namespace Utils