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

#ifndef APPLICATION_LAYERIO_ITKLAYERIMPORTER_H
#define APPLICATION_LAYERIO_ITKLAYERIMPORTER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

// Application includes
#include <Application/LayerIO/LayerImporter.h>
#include <Application/LayerIO/LayerIO.h>

namespace Seg3D
{

class ITKLayerImporter : public LayerImporter
{
	SCI_IMPORTER_TYPE( "ITK Importer", "*", 5)

	// -- Constructor/Destructor --
public:
	// Construct a new layer file importer
	ITKLayerImporter(const std::string& filename) :
		LayerImporter(filename)
	{
	}

	// Virtual destructor for memory management of derived classes
	virtual ~ITKLayerImporter()
	{
	}

	// -- Import a file information --
public:

	// IMPORT_HEADER:
	// Import all the information needed to translate the header and metadata information, but not
	// necessarily read the whole file. NOTE: Some external packages do not support reading a header
	// and hence these importers should read the full file here.
	virtual bool import_header();

	// GET_GRID_TRANSFORM:
	// Get the grid transform of the grid that we are importing
	virtual Utils::GridTransform get_grid_transform();

	// GET_DATA_TYPE:
	// Get the type of data that is being imported
	virtual Utils::DataType get_data_type();

	// GET_IMPORTER_MODES:
	// Get then supported importer modes
	virtual int get_importer_modes();
	
	// --Import the data as a specific type --	
public:	

	// IMPORT_LAYER
	// Import the layer from the file
	virtual bool import_layer( LayerImporterMode mode, std::vector<LayerHandle>& layers );

};

} // end namespace seg3D

#endif
