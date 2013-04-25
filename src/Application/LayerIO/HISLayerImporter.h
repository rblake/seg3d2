/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#ifndef APPLICATION_LAYERIO_HISLAYERIMPORTER_H
#define APPLICATION_LAYERIO_HISLAYERIMPORTER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// Application includes
#include <Application/LayerIO/LayerFileSeriesImporter.h>
#include <Application/LayerIO/LayerIO.h>

namespace Seg3D
{

// Forward declaration for internal class
class HISLayerImporterPrivate;
typedef boost::shared_ptr<HISLayerImporterPrivate> HISLayerImporterPrivateHandle;  

class HISLayerImporter : public LayerFileSeriesImporter
{
  SEG3D_IMPORTER_TYPE( "HIS Importer", ".his;.HIS", 80 )

public:
  HISLayerImporter();
  virtual ~HISLayerImporter();
  
  // -- Import information from file --
public:
  /// GET_FILE_INFO
  /// Get the information about the file we are currently importing.
  /// NOTE: This function often causes the file to be loaded in its entirety
  /// Hence it is best to run this on a separate thread if needed ( from the GUI ).
  virtual bool get_file_info( LayerImporterFileInfoHandle& info );
  
  // -- Import data from file --	
public:	
  /// GET_FILE_DATA
  /// Get the file data from the file/ file series
  /// NOTE: The information is generated again, so that hints can be processed
  virtual bool get_file_data( LayerImporterFileDataHandle& data );
  
  // --internals --
public:
  HISLayerImporterPrivateHandle private_;
};

}

#endif