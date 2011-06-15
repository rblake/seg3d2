/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2011 Scientific Computing and Imaging Institute,
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
 
/*
 * MRC file format reader: http://www2.mrc-lmb.cam.ac.uk/image2000.html
 * Implementation follows EMAN2 and Chimera.
 */


#include "MRCReader.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace MRC2000IO {

bool isBigEndian()
{
  unsigned short i = 0x4321;
  if ( (*(reinterpret_cast<unsigned char *>(&i))) != 0x21 )
  {
    return true;
  } 
  else
  {
    return false;
  }
}
  
MRCReader::MRCReader()
  : host_is_big_endian_(isBigEndian()),
    swap_endian_(false),
    use_new_origin_(true),
    MASK_UPPER_(0xFF000000),
    MASK_LOWER_(0x000000FF)
{
}

MRCReader::~MRCReader() {}


bool MRCReader::read_header(const std::string& filename, MRCHeader& header)
{
  try
  {
    header.machinestamp = 0;
    header.map[3] = '\0';

    std::ifstream::pos_type size = 0;

    std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
    if (! in)
    {
      error_ = std::string("Failed to open file ") + filename;
      return false;
    }

    in.seekg(0, std::ios::end);
    size = in.tellg();
    std::cout << "File is " << size << " bytes long." << std::endl;
    in.seekg(0, std::ios::beg);
    std::ifstream::pos_type data_block_size = size - static_cast<std::ifstream::pos_type>(MRC_HEADER_LENGTH);
    std::cout << "Data block is " << data_block_size << " bytes long." << std::endl;

    in.read(reinterpret_cast<char*>(&header), MRC_HEADER_LENGTH);
    process_header(&header, MRC_HEADER_LENGTH);

    if (header.mode == MRC_SHORT_COMPLEX || header.mode == MRC_FLOAT_COMPLEX)
    {
      error_ = "Complex mode is not supported.";
      return false;
    }

    // test
    std::cout << "new origin: [" << header.xorigin << " " << header.yorigin << " " << header.zorigin << "]" << std::endl;
    std::cout << "old origin: [" << header.nxstart << " " << header.nystart << " " << header.nzstart << "]" << std::endl;
    std::cout << "cell dims in angstroms: [" << header.xlen << " " << header.ylen << " " << header.zlen << "]" << std::endl;
    // test
    if ((0 == header.xorigin || 0 == header.yorigin || 0 == header.zorigin) &&
        (0 != header.nxstart || 0 != header.nystart || 0 != header.nzstart))
    {
      // use n[x|y|z]start
      this->use_new_origin_ = false;
    }
    in.close();

    // test
    std::hex(std::cout);
    std::cout << "header machinestamp=" << header.machinestamp << std::endl;
    std::dec(std::cout);
    std::cout << header.map << std::endl;
    // test
  }
  // TODO: ios specific exceptions
  catch (...)
  {
    error_ = std::string("Failed to read header for file ") + filename;
    return false;
  }

  return true;
}

bool MRCReader::process_header(void* buffer, int buffer_len)
{
  int* long_word_buffer = reinterpret_cast<int*>(buffer);
  const int machine_stamp = long_word_buffer[53];

  // file endianness vs. host endianness
  // N.B. machine stamp field is not always implemented reliably
  if (host_is_big_endian_)
  {
    // big endian host, data read on little endian machine
    if (( (machine_stamp & MASK_LOWER_) == 0x44 ) || ( (machine_stamp & MASK_UPPER_) == 0x44000000 ))
    {
      swap_endian_ = true;
    }
    // big endian host, data read on big endian machine
    else if ( (machine_stamp &  MASK_LOWER_) == 0x11)
    {
      swap_endian_ = false;
    }
    else
    {
      // Guess endianness of data by checking upper bytes of nx:
      // assumes nx < 2^16
      char* char_buffer = reinterpret_cast<char*>(buffer);
      swap_endian_ = true;
      for (size_t j = MRC_LONG_WORD / 2; j < MRC_LONG_WORD; ++j) {
        // test
        //std::cout << static_cast<int>(char_buffer[j]) << std::endl;
        // test
        if (char_buffer[j] != 0) {
          swap_endian_ = false;
          break;
        }
      }
    }
  }
  else // little endian (no other architectures supported by Core/Endian)
  {
    // little endian host, data read on little endian machine
    if (( (machine_stamp & MASK_LOWER_) == 0x44 ) || ( (machine_stamp & MASK_UPPER_) == 0x44000000 ))
    {
      swap_endian_ = false;
    }
    // little endian host, data read on big endian machine
    else if ( (machine_stamp &  MASK_LOWER_) == 0x11)
    {
      swap_endian_ = true;
    }
    else
    {
      // Guess endianness of data by checking lower bytes of nx:
      // assumes nx < 2^16
      char* char_buffer = reinterpret_cast<char*>(buffer);
      swap_endian_ = true;
      for (size_t j = 0; j < MRC_LONG_WORD / 2; ++j) {
        // test
        //std::cout << static_cast<int>(char_buffer[j]) << std::endl;
        // test
        if (char_buffer[j] != 0) {
          swap_endian_ = false;
          break;
        }
      }
    }
  }

  // convert header to little endian before reading
  if (swap_endian_)
  {
    // code from DataBlock.cc, SwapEndian method
    unsigned char* ubuffer = reinterpret_cast<unsigned char*>(long_word_buffer);
    unsigned char tmp;
    const size_t SIZE8 = MRC_HEADER_LENGTH_LWORDS & ~(0x7);
    // test
    //std::hex(std::cout);
    //std::cout << "SIZE8=" << SIZE8 << std::endl;
    // test
    size_t i = 0;
    // swap word bytes in blocks of 32 bytes in place
    for(; i < SIZE8; ++i)
    {
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;	

      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
    }
    //std::dec(std::cout);
    //std::cout << "i=" << i << std::endl;

    for(; i < MRC_HEADER_LENGTH_LWORDS; ++i)
    {
      tmp = ubuffer[ 0 ]; ubuffer[ 0 ] = ubuffer[ 3 ]; ubuffer[ 3 ] = tmp;
      tmp = ubuffer[ 1 ]; ubuffer[ 1 ] = ubuffer[ 2 ]; ubuffer[ 2 ] = tmp;
      ubuffer += 4;
    }
  }

  MRCHeader *h = reinterpret_cast<MRCHeader*>(long_word_buffer);
  int temp = 0;

  temp = long_word_buffer[10];
  h->xlen = static_cast<float>(temp);
  temp = long_word_buffer[11];
  h->ylen = static_cast<float>(temp);
  temp = long_word_buffer[12];
  h->zlen = static_cast<float>(temp);

  temp = long_word_buffer[13];
  h->alpha = static_cast<float>(temp);
  temp = long_word_buffer[14];
  h->beta = static_cast<float>(temp);
  temp = long_word_buffer[15];
  h->gamma = static_cast<float>(temp);

  temp = long_word_buffer[19];
  h->dmin = static_cast<float>(temp);
  temp = long_word_buffer[20];
  h->dmax = static_cast<float>(temp);
  temp = long_word_buffer[21];
  h->dmean = static_cast<float>(temp);
  temp = long_word_buffer[22];
  h->ispg = static_cast<int>(temp);
  temp = long_word_buffer[23];
  h->nsymbt = static_cast<int>(temp);

  temp = long_word_buffer[49];
  h->xorigin = static_cast<float>(temp);
  temp = long_word_buffer[50];
  h->yorigin = static_cast<float>(temp);
  temp = long_word_buffer[51];
  h->zorigin = static_cast<float>(temp);
  // Force a sentinal null char
  h->map[MRC_LONG_WORD-1] = '\0';

  temp = long_word_buffer[54];
  h->rms = static_cast<float>(temp);

  return true;
}

}
