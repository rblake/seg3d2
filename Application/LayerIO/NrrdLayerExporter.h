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

#ifndef APPLICATION_LAYERIO_NRRDLAYEREXPORTER_H
#define APPLICATION_LAYERIO_NRRDLAYEREXPORTER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 

#include <boost/filesystem.hpp>


// Application includes
#include <Application/LayerIO/LayerExporter.h>
#include <Application/LayerIO/LayerIO.h>

namespace Seg3D
{

class NrrdLayerExporter : public LayerExporter
{
	SCI_EXPORTER_TYPE("NRRD Exporter", ".nrrd")

	// -- Constructor/Destructor --
public:
	// Construct a new layer file importer
	NrrdLayerExporter( std::vector< LayerHandle >& layers );

	// Virtual destructor for memory management of derived classes
	virtual ~NrrdLayerExporter()
	{
	}

	// -- Import a file information --
public:
	// GET_GRID_TRANSFORM:
	// Get the grid transform of the grid that we are importing
	virtual Core::GridTransform get_grid_transform();

	// GET_DATA_TYPE:
	// Get the type of data that is being imported
	virtual Core::DataType get_data_type();

	// GET_IMPORTER_MODES
	// Get then supported importer modes
	virtual int get_exporter_modes();
	
	// --Import the data as a specific type --	
public:	

	// EXPORT_LAYER
	// Export the layer to file
	virtual bool export_layer( LayerExporterMode mode, const std::string& file_path, 
		const std::string& name );
		
	virtual void set_label_layer_values( std::vector< double > values ){ this->label_values_ = values; }
		
private:
	bool export_nrrd( const std::string& file_path );
	bool export_single_masks( const std::string& file_path );
	bool export_mask_label( const std::string& file_path );
	
private:
	std::vector< double > label_values_;

};

} // end namespace seg3D

#endif
