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

// Application includes
#include <Application/LayerIO/ITKLayerImporter.h>

SCI_REGISTER_IMPORTER( Seg3D, ITKLayerImporter );


namespace Seg3D
{

bool ITKLayerImporter::import_header()
{
	return false;
}

Core::GridTransform ITKLayerImporter::get_grid_transform()
{
	return Core::GridTransform(1,1,1);
}

Core::DataType ITKLayerImporter::get_data_type()
{
	return Core::DataType::UNKNOWN_E;
}

int ITKLayerImporter::get_importer_modes()
{
	return 0;
}

bool ITKLayerImporter::import_layer( LayerImporterMode mode, std::vector<LayerHandle>& layers )
{
	return false;
}

} // end namespace seg3D
