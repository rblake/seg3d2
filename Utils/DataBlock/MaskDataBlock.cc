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

#include <Utils/DataBlock/MaskDataBlock.h>
#include <Utils/DataBlock/MaskDataBlockManager.h>

namespace Utils
{

MaskDataBlock::MaskDataBlock( DataBlockHandle& data_block, unsigned int mask_bit ) :
	nx_( data_block->get_nx() ),
	ny_( data_block->get_ny() ),
	nz_( data_block->get_nz() ),
	data_block_( data_block ),
	mask_bit_( mask_bit ),
	bit_tester_( 1 << mask_bit )
{
	this->data_ = reinterpret_cast<unsigned char*>( this->data_block_->get_data() );
}

MaskDataBlock::~MaskDataBlock()
{
	MaskDataBlockManager::Instance()->release( data_block_, mask_bit_ );
}

} // end namespace Utils
