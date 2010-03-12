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
#include <Application/LayerIO/NrrdLayerImporter.h>

namespace Seg3D
{

SCI_REGISTER_IMPORTER(NrrdLayerImporter);


bool NrrdLayerImporter::import_header()
{
	return false;
}

bool NrrdLayerImporter::import_data()
{
	return false;
}

Utils::GridTransform NrrdLayerImporter::get_grid_transform()
{
	Utils::GridTransform identity(1,1,1);
	return identity;
}

bool NrrdLayerImporter::is_data_volume_compatible()
{
	return false;
}

bool NrrdLayerImporter::is_mask_volume_compatible()
{
	return false;
}

bool NrrdLayerImporter::is_label_volume_compatible()
{
	return false;
}

bool NrrdLayerImporter::import_as_datavolume( LayerHandle& layer )
{
	return false;
}

bool NrrdLayerImporter::import_as_maskvolume( std::vector<LayerHandle>& layers,
		LayerMaskImporterMode mode )
{
	return false;
}

bool NrrdLayerImporter::import_as_labelvolume( LayerHandle& layer )
{
	return false;
}

} // end namespace seg3D
