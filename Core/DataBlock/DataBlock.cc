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

#include <Core/DataBlock/DataBlock.h>
#include <Core/DataBlock/DataBlockManager.h>
#include <Core/DataBlock/StdDataBlock.h>

namespace Core
{

DataBlock::DataBlock() :
	nx_( 0 ), 
	ny_( 0 ), 
	nz_( 0 ), 
	data_type_( DataType::UNKNOWN_E ), 
	data_( 0 ),
	generation_( -1 )
{
}

DataBlock::~DataBlock()
{
	// Unregister the data block if it's registered
	if ( this->generation_ != -1 )
	{
		DataBlockManager::Instance()->unregister_datablock( this->generation_ );
	}
}

double DataBlock::get_data_at( size_t index ) const
{
	switch( this->data_type_ )
	{
	case DataType::CHAR_E:
		{
			signed char* data = reinterpret_cast<signed char*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}			
	case DataType::UCHAR_E:
		{
			unsigned char* data = reinterpret_cast<unsigned char*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}			
	case DataType::SHORT_E:
		{
			short* data = reinterpret_cast<short*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}			
	case DataType::USHORT_E:
		{
			unsigned short* data = reinterpret_cast<unsigned short*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}
	case DataType::INT_E:
		{
			int* data = reinterpret_cast<int*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}			
	case DataType::UINT_E:
		{
			unsigned int* data = reinterpret_cast<unsigned int*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}				
	case DataType::FLOAT_E:
		{
			float* data = reinterpret_cast<float*>( this->data_ );
			return static_cast<double>( data[ index ] );
		}			
	case DataType::DOUBLE_E:
		{
			double* data = reinterpret_cast<double*>( this->data_ );
			return data[ index ];
		}			
	}
	
	return 0.0;
}

void DataBlock::set_data_at( size_t index, double value )
{
	switch( this->data_type_ )
	{
		case DataType::CHAR_E:
		{
			signed char* data = reinterpret_cast<signed char*>( this->data_ );
			data[ index ] = static_cast<signed char>( value );
			return;
		}			
		case DataType::UCHAR_E:
		{
			unsigned char* data = reinterpret_cast<unsigned char*>( this->data_ );
			data[ index ] = static_cast<unsigned char>( value );
			return;
		}			
		case DataType::SHORT_E:
		{
			short* data = reinterpret_cast<short*>( this->data_ );
			data[ index ] = static_cast<short>( value );
			return;
		}			
		case DataType::USHORT_E:
		{
			unsigned short* data = reinterpret_cast<unsigned short*>( this->data_ );
			data[ index ] = static_cast<unsigned short>( value );
			return;
		}	
		case DataType::INT_E:
		{
			int* data = reinterpret_cast<int*>( this->data_ );
			data[ index ] = static_cast<int>( value );
			return;
		}			
		case DataType::UINT_E:
		{
			unsigned int* data = reinterpret_cast<unsigned int*>( this->data_ );
			data[ index ] = static_cast<unsigned int>( value );
			return;
		}	
		case DataType::FLOAT_E:
		{
			float* data = reinterpret_cast<float*>( this->data_ );
			data[ index ] = static_cast<float>( value );
			return;
		}			
		case DataType::DOUBLE_E:
		{
			double* data = reinterpret_cast<double*>( this->data_ );
			data[ index ] = value;
			return;
		}
	}
}

void DataBlock::set_type( DataType type )
{
	this->data_type_ = type;
}

void DataBlock::set_nx( size_t nx )
{
	this->nx_ = nx;
}

void DataBlock::set_ny( size_t ny )
{
	this->ny_ = ny;
}

void DataBlock::set_nz( size_t nz )
{
	this->nz_ = nz;
}

void DataBlock::set_data( void* data )
{
	this->data_ = data;
}

void DataBlock::set_histogram( const Histogram& histogram )
{
	this->histogram_ = histogram;
}

double DataBlock::get_max() const
{
	return this->histogram_.get_max();
}

double DataBlock::get_min() const
{
	return this->histogram_.get_min();
}

double DataBlock::get_range() const
{
	return this->histogram_.get_max() - this->histogram_.get_min();
}

const Histogram& DataBlock::get_histogram() const
{
	return this->histogram_;
}

void DataBlock::clear()
{
	lock_type lock( this->get_mutex() );
	memset( this->data_, 0, Core::GetSizeDataType( this->data_type_ ) * this->get_size() );
	this->generation_ = DataBlockManager::Instance()->increase_generation( this->generation_ );
}

DataBlock::generation_type DataBlock::get_generation() const
{
	shared_lock_type lock( this->get_mutex() );
	return this->generation_;
}

void DataBlock::set_generation( generation_type generation )
{
	lock_type lock( this->get_mutex() );
	this->generation_ = generation;
}

void DataBlock::increase_generation()
{
	lock_type lock( this->get_mutex() );
	this->generation_ = DataBlockManager::Instance()->increase_generation( this->generation_ );
}

bool DataBlock::update_histogram()
{
	lock_type lock( this->get_mutex() );

	switch( this->data_type_ )
	{
		case DataType::CHAR_E:
			return this->histogram_.compute( reinterpret_cast<signed char*>( get_data() ), get_size() );
		case DataType::UCHAR_E:
			return this->histogram_.compute( reinterpret_cast<unsigned char*>( get_data() ), get_size() );
		case DataType::SHORT_E:
			return this->histogram_.compute( reinterpret_cast<short*>( get_data() ), get_size() );
		case DataType::USHORT_E:
			return this->histogram_.compute( reinterpret_cast<unsigned short*>( get_data() ), get_size() );
		case DataType::INT_E:
			return this->histogram_.compute( reinterpret_cast<int*>( get_data() ), get_size() );
		case DataType::UINT_E:
			return this->histogram_.compute( reinterpret_cast<unsigned int*>( get_data() ), get_size() );
		case DataType::FLOAT_E:
			return this->histogram_.compute( reinterpret_cast<float*>( get_data() ), get_size() );
		case DataType::DOUBLE_E:
			return this->histogram_.compute( reinterpret_cast<double*>( get_data() ), get_size() );
	}

	return false;
}

void DataBlock::update_data_type( DataType type )
{
	this->data_type_ = type;
}

template<class DATA>
static bool ConvertDataTypeInternal( DATA* src, DataBlockHandle& dst_data_block )
{
	size_t size = dst_data_block->get_size();
	switch ( dst_data_block->get_data_type() )
	{
		case DataType::CHAR_E:	
		{
			signed char* dst = reinterpret_cast<signed char*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<signed char>( src[ j ] );
				dst[ j + 1 ] = static_cast<signed char>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<signed char>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<signed char>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<signed char>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<signed char>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<signed char>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<signed char>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<signed char>( src[ j ] );
			}
			return true;
		}
		case DataType::UCHAR_E:	
		{
			unsigned char* dst = reinterpret_cast<unsigned char*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<unsigned char>( src[ j ] );
				dst[ j + 1 ] = static_cast<unsigned char>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<unsigned char>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<unsigned char>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<unsigned char>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<unsigned char>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<unsigned char>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<unsigned char>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<unsigned char>( src[ j ] );
			}			
			return true;
		}
		case DataType::SHORT_E:	
		{
			short* dst = reinterpret_cast<short*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<short>( src[ j ] );
				dst[ j + 1 ] = static_cast<short>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<short>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<short>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<short>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<short>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<short>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<short>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<short>( src[ j ] );
			}			
			return true;
		}
		case DataType::USHORT_E:	
		{
			unsigned short* dst = reinterpret_cast<unsigned short*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<unsigned short>( src[ j ] );
				dst[ j + 1 ] = static_cast<unsigned short>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<unsigned short>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<unsigned short>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<unsigned short>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<unsigned short>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<unsigned short>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<unsigned short>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<unsigned short>( src[ j ] );
			}	
			return true;
		}
		case DataType::INT_E:	
		{
			int* dst = reinterpret_cast<int*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<int>( src[ j ] );
				dst[ j + 1 ] = static_cast<int>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<int>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<int>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<int>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<int>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<int>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<int>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<int>( src[ j ] );
			}		
			return true;
		}
		case DataType::UINT_E:	
		{
			unsigned int* dst = reinterpret_cast<unsigned int*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<unsigned int>( src[ j ] );
				dst[ j + 1 ] = static_cast<unsigned int>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<unsigned int>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<unsigned int>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<unsigned int>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<unsigned int>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<unsigned int>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<unsigned int>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<unsigned int>( src[ j ] );
			}	
			return true;
		}
		case DataType::FLOAT_E:	
		{
			float* dst = reinterpret_cast<float*>( dst_data_block->get_data() );
			size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<float>( src[ j ] );
				dst[ j + 1 ] = static_cast<float>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<float>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<float>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<float>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<float>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<float>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<float>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<float>( src[ j ] );
			}				
			return true;
		}
		case DataType::DOUBLE_E:	
		{
			double* dst = reinterpret_cast<double*>( dst_data_block->get_data() );
						size_t size8 = size & ~(0x7);
			size_t j = 0;
			for ( ; j < size8; j+=8 )
			{
				dst[ j ] = static_cast<double>( src[ j ] );
				dst[ j + 1 ] = static_cast<double>( src[ j + 1 ] );
				dst[ j + 2 ] = static_cast<double>( src[ j + 2 ] );
				dst[ j + 3 ] = static_cast<double>( src[ j + 3 ] );
				dst[ j + 4 ] = static_cast<double>( src[ j + 4 ] );
				dst[ j + 5 ] = static_cast<double>( src[ j + 5 ] );
				dst[ j + 6 ] = static_cast<double>( src[ j + 6 ] );
				dst[ j + 7 ] = static_cast<double>( src[ j + 7 ] );
			}
			for ( ; j < size; j++ )
			{
				dst[ j ] = static_cast<double>( src[ j ] );
			}	
			return true;
		}
		default:
		{
			dst_data_block.reset();
			return false;
		}
	}
}


bool DataBlock::ConvertDataType( const DataBlockHandle& src_data_block, 
	DataBlockHandle& dst_data_block, DataType new_data_type )
{
	dst_data_block.reset();
	if ( !src_data_block ) return false;

	shared_lock_type lock( src_data_block->get_mutex( ) );

	dst_data_block = StdDataBlock::New( src_data_block->get_nx(),
		src_data_block->get_ny(), src_data_block->get_nz(), new_data_type );
		
	if ( !dst_data_block )
	{
		return false;
	}
	
	switch( src_data_block->get_data_type() )
	{
		case DataType::CHAR_E:
			return ConvertDataTypeInternal<signed char>( 
				reinterpret_cast<signed char*>( src_data_block->get_data() ), dst_data_block );
		case DataType::UCHAR_E:
			return ConvertDataTypeInternal<unsigned char>( 
				reinterpret_cast<unsigned char*>( src_data_block->get_data() ), dst_data_block );
		case DataType::SHORT_E:
			return ConvertDataTypeInternal<short>( 
				reinterpret_cast<short*>( src_data_block->get_data() ), dst_data_block );
		case DataType::USHORT_E:
			return ConvertDataTypeInternal<unsigned short>( 
				reinterpret_cast<unsigned short*>( src_data_block->get_data() ), dst_data_block );
		case DataType::INT_E:
			return ConvertDataTypeInternal<int>( 
				reinterpret_cast<int*>( src_data_block->get_data() ), dst_data_block );
		case DataType::UINT_E:
			return ConvertDataTypeInternal<unsigned int>( 
				reinterpret_cast<unsigned int*>( src_data_block->get_data() ), dst_data_block );
		case DataType::FLOAT_E:
			return ConvertDataTypeInternal<float>( 
				reinterpret_cast<float*>( src_data_block->get_data() ), dst_data_block );
		case DataType::DOUBLE_E:
			return ConvertDataTypeInternal<double>( 
				reinterpret_cast<double*>( src_data_block->get_data() ), dst_data_block );
		default:
			dst_data_block.reset();
			return false;
	}
}

template<class DATA>
static bool PermuteDataInternal( const DataBlockHandle& src_data_block, 
	DataBlockHandle& dst_data_block, std::vector<int>& permutation )
{
	DATA* src = reinterpret_cast<DATA*>( src_data_block->get_data() );
	DATA* dst = reinterpret_cast<DATA*>( dst_data_block->get_data() );

	typedef DataBlock::index_type index_type;

	std::vector<index_type> start(3);
	std::vector<index_type> stride(3);
	
	index_type nx = static_cast<index_type>( src_data_block->get_nx() );
	index_type ny = static_cast<index_type>( src_data_block->get_ny() );
	index_type nz = static_cast<index_type>( src_data_block->get_nz() );
	index_type nxy = nx*ny;
	index_type nxyz = nx*ny*nz;
	
	for ( index_type j = 0; j < 3; j++)
	{
		if ( permutation[ j ] == 1 )
		{
			start[ j ] = 0;
			stride[ j ] = 1; 
		}
		else if ( permutation[ j ] == -1 )
		{
			start[ j ] = nx - 1;
			stride[ j ] = -1; 	
		}
		else if ( permutation[ j ] == 2 )
		{
			start[ j ] = 0;
			stride[ j ] = nx; 
		}
		else if ( permutation[ j ] == -2 )
		{
			start[ j ] = nxy - nx;
			stride[ j ] = -nx; 	
		}
		else if ( permutation[ j ] == 3 )
		{
			start[ j ] = 0;
			stride[ j ] = nxy; 
		}
		else if ( permutation[ j ] == -3 )
		{
			start[ j ] = nxyz - nxy;
			stride[ j ] = -nxy; 	
		}
	}
	
	index_type dnx = static_cast<index_type>( dst_data_block->get_nx() );
	index_type dny = static_cast<index_type>( dst_data_block->get_ny() );
	index_type dnz = static_cast<index_type>( dst_data_block->get_nz() );
	index_type dnxy = dnx * dny;
	index_type dnxyz = dnx * dny * dnz;
	index_type sz = start[ 2 ];
	index_type sz_stride = stride[ 2 ];
	for ( index_type dz = 0; dz < dnxyz; dz += dnxy )
	{
		index_type sy = start[ 1 ];
		index_type sy_stride = stride[ 1 ];
		for ( index_type dy = 0; dy < dnxy; dy += dnx )
		{
			index_type sx = start[ 0 ];
			index_type sx_stride = stride[ 0 ];
			for ( index_type dx = 0; dx < dnx; dx++ )
			{
				dst[ dx + dy + dz ] = src[ sx + sy + sz ];
				sx += sx_stride;
			}
			sy += sy_stride;
		}
		sz += sz_stride;
	}
	
	return true;
}

bool DataBlock::PermuteData( const DataBlockHandle& src_data_block, 
	DataBlockHandle& dst_data_block, std::vector<int> permutation )
{
	dst_data_block.reset();
	if ( !src_data_block ) return false;

	shared_lock_type lock( src_data_block->get_mutex( ) );

	if ( permutation.size() != 3 )
	{
		dst_data_block.reset();
		return false;
	}
	
	size_t dn[3];
	dn[ 0 ] = 0;
	dn[ 1 ] = 0;
	dn[ 2 ] = 0;
	
	for ( size_t j = 0; j < 3; j++)
	{
		if ( permutation[ j ] == 1 || permutation[ j ] == -1 ) dn[ j ] = src_data_block->get_nx();
		if ( permutation[ j ] == 2 || permutation[ j ] == -2 ) dn[ j ] = src_data_block->get_ny();
		if ( permutation[ j ] == 3 || permutation[ j ] == -3 ) dn[ j ] = src_data_block->get_nz();		
	}
	
	if ( dn[ 0 ] * dn[ 1 ] * dn[ 2 ] == 0 ) return false;
	
	dst_data_block = StdDataBlock::New( dn[ 0 ], dn[ 1 ], dn[ 2 ], 
		src_data_block->get_data_type() );	

	switch( src_data_block->get_data_type() )
	{
		case DataType::CHAR_E:
			return PermuteDataInternal<signed char>( src_data_block, dst_data_block, 
				permutation );
		case DataType::UCHAR_E:
			return PermuteDataInternal<unsigned char>( src_data_block, dst_data_block, 
				permutation );
		case DataType::SHORT_E:
			return PermuteDataInternal<short>( src_data_block, dst_data_block, 
				permutation );
		case DataType::USHORT_E:
			return PermuteDataInternal<unsigned short>( src_data_block, dst_data_block, 
				permutation );
		case DataType::INT_E:
			return PermuteDataInternal<int>( src_data_block, dst_data_block, 
				permutation );
		case DataType::UINT_E:
			return PermuteDataInternal<unsigned int>( src_data_block, dst_data_block, 
				permutation );
		case DataType::FLOAT_E:
			return PermuteDataInternal<float>( src_data_block, dst_data_block, 
				permutation );
		case DataType::DOUBLE_E:
			return PermuteDataInternal<double>( src_data_block, dst_data_block, 
				permutation );
		default:
			return false;
	}
}


template<class DATA>
static bool QuantizeDataInternal( double min, double max, DATA* src, DataBlockHandle& dst_data_block )
{
	float fmin = static_cast<float>( min );
	float fmax = static_cast<float>( max );
	
	size_t size = dst_data_block->get_size();
	switch ( dst_data_block->get_data_type() )
	{
		case DataType::CHAR_E:	
		{
			signed char* dst = reinterpret_cast<signed char*>( dst_data_block->get_data() );

			float offset = 0.5f - static_cast<float>( 0x80 );
			float multiplier = 0.0f;
			if ( fmax > fmin ) multiplier = static_cast<float>( 0x100 ) / (fmax - fmin);
			
			for ( size_t j = 0 ; j < size; j++ )
			{
				dst[ j ] = static_cast<signed char>( multiplier * 
					(static_cast<float>( src[ j ] ) - fmin)  + offset );
			}
			
			return true;
		}
		case DataType::UCHAR_E:	
		{
			unsigned char* dst = reinterpret_cast<unsigned char*>( dst_data_block->get_data() );
			
			float offset = 0.5f;
			float multiplier = 0.0f;
			if ( fmax > fmin ) multiplier = static_cast<float>( 0x100 ) / (fmax - fmin);
			 			 
			for ( size_t j = 0 ; j < size; j++ )
			{
				dst[ j ] = static_cast<unsigned char>( multiplier * 
					(static_cast<float>( src[ j ] ) - fmin)  + offset );
			}
			return true;
		}
		case DataType::SHORT_E:	
		{
			short* dst = reinterpret_cast<short*>( dst_data_block->get_data() );

			float offset = 0.5f - static_cast<float>( 0x8000 );
			float multiplier = 0.0f;
			if ( fmax > fmin ) multiplier = static_cast<float>( 0x10000 ) / (fmax - fmin);
			
			for ( size_t j = 0 ; j < size; j++ )
			{
				dst[ j ] = static_cast<short>( multiplier * 
					(static_cast<float>( src[ j ] ) - fmin)  + offset );
			}
			
			return true;
		}
		case DataType::USHORT_E:	
		{
			unsigned short* dst = reinterpret_cast<unsigned short*>( dst_data_block->get_data() );

			float offset = 0.5f;
			float multiplier = 0.0f;
			if ( fmax > fmin ) multiplier = static_cast<float>( 0x10000 ) / (fmax - fmin);
			
			for ( size_t j = 0 ; j < size; j++ )
			{
				dst[ j ] = static_cast<unsigned short>( multiplier * 
					(static_cast<float>( src[ j ] ) - fmin)  + offset);
			}
			
			return true;
		}
		case DataType::INT_E:	
		{
			int* dst = reinterpret_cast<int*>( dst_data_block->get_data() );

			double offset = 0.5 -  static_cast<double>( 0x80000000 );
			double multiplier = 0.0;
			if ( max > min ) multiplier = static_cast<double>( 0x100000000ull ) / (max - min);
			
			for ( size_t j = 0 ; j < size; j++ )
			{
				dst[ j ] = static_cast<int>( multiplier * 
					(static_cast<double>( src[ j ] ) - min)  + offset);
			}
			
			return true;
		}
		case DataType::UINT_E:	
		{
			unsigned int* dst = reinterpret_cast<unsigned int*>( dst_data_block->get_data() );

			double offset = 0.5;
			double multiplier = 0.0;
			if ( max > min ) multiplier = static_cast<double>( 0x100000000ull ) / (max - min);
			
			for ( size_t j = 0 ; j < size; j++ )
			{
				dst[ j ] = static_cast<unsigned int>( multiplier * 
					(static_cast<double>( src[ j ] ) - min)  + offset);
			}
			
			return true;
		}
		default:
		{
			dst_data_block.reset();
			return false;
		}
	}
}

bool DataBlock::QuantizeData( const DataBlockHandle& src_data_block, 
	DataBlockHandle& dst_data_block, DataType new_data_type )
{
	dst_data_block.reset();
	if ( !src_data_block ) return false;
	
	shared_lock_type lock( src_data_block->get_mutex( ) );

	if ( new_data_type != DataType::CHAR_E && new_data_type != DataType::UCHAR_E &&
		new_data_type != DataType::SHORT_E && new_data_type != DataType::USHORT_E &&
		new_data_type != DataType::INT_E && new_data_type != DataType::UINT_E )
	{
		return false;
	}	

	dst_data_block = StdDataBlock::New( src_data_block->get_nx(),
		src_data_block->get_ny(), src_data_block->get_nz(), new_data_type );
		
	if ( !dst_data_block )
	{
		return false;
	}
	
	double min = src_data_block->get_min();
	double max = src_data_block->get_max();
	
	switch( src_data_block->get_data_type() )
	{
		case DataType::CHAR_E:
			return QuantizeDataInternal<signed char>( min, max, 
				reinterpret_cast<signed char*>( src_data_block->get_data() ), dst_data_block );
		case DataType::UCHAR_E:
			return QuantizeDataInternal<unsigned char>( min, max, 
				reinterpret_cast<unsigned char*>( src_data_block->get_data() ), dst_data_block );
		case DataType::SHORT_E:
			return QuantizeDataInternal<short>( min, max,
				reinterpret_cast<short*>( src_data_block->get_data() ), dst_data_block );
		case DataType::USHORT_E:
			return QuantizeDataInternal<unsigned short>( min, max,
				reinterpret_cast<unsigned short*>( src_data_block->get_data() ), dst_data_block );
		case DataType::INT_E:
			return QuantizeDataInternal<int>( min, max,
				reinterpret_cast<int*>( src_data_block->get_data() ), dst_data_block );
		case DataType::UINT_E:
			return QuantizeDataInternal<unsigned int>( min, max,
				reinterpret_cast<unsigned int*>( src_data_block->get_data() ), dst_data_block );
		case DataType::FLOAT_E:
			return QuantizeDataInternal<float>( min, max,
				reinterpret_cast<float*>( src_data_block->get_data() ), dst_data_block );
		case DataType::DOUBLE_E:
			return QuantizeDataInternal<double>( min, max,
				reinterpret_cast<double*>( src_data_block->get_data() ), dst_data_block );
		default:
			return false;
	}
}

bool DataBlock::Clone( const DataBlockHandle& src_data_block, 
		DataBlockHandle& dst_data_block )
{
	// Step (1) : Check whether there is a source data block
	dst_data_block.reset();
	if ( !src_data_block ) return false;

	// Step (2) : Lock the source
	shared_lock_type lock( src_data_block->get_mutex( ) );

	// Step (3): Generate a new data block with the right type
	dst_data_block = StdDataBlock::New( src_data_block->get_nx(),
		src_data_block->get_ny(), src_data_block->get_nz(), src_data_block->get_data_type() );
		
	// Step (4): Copy the data	
	size_t mem_size = src_data_block->get_size();	
	switch( src_data_block->get_data_type() )
	{
		case DataType::CHAR_E:
			mem_size *= sizeof( signed char );
			break;
		case DataType::UCHAR_E:
			mem_size *= sizeof( unsigned char );
			break;
		case DataType::SHORT_E:
			mem_size *= sizeof( short );
			break;
		case DataType::USHORT_E:
			mem_size *= sizeof( unsigned short );
			break;
		case DataType::INT_E:
			mem_size *= sizeof( int );
			break;
		case DataType::UINT_E:
			mem_size *= sizeof( unsigned int );
			break;
		case DataType::FLOAT_E:
			mem_size *= sizeof( float );
			break;
		case DataType::DOUBLE_E:
			mem_size *= sizeof( double );
			break;
		default:
			return false;
	}
	std::memcpy( dst_data_block->get_data(), src_data_block->get_data(), mem_size );
	
	// Step (5) : Copy the histogram
	dst_data_block->set_histogram( src_data_block->get_histogram() );

	return true;
}

} // end namespace Core
