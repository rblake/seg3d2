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

#ifndef APPLICATION_LAYERMANAGER_LAYERSCENE_H
#define APPLICATION_LAYERMANAGER_LAYERSCENE_H

#include <vector>
#include <boost/shared_ptr.hpp>

#include <Application/Layer/Layer.h>

#include <Utils/Volume/DataVolume.h>
#include <Utils/Volume/MaskVolume.h>

namespace Seg3D
{

class LayerSceneItem;
typedef boost::shared_ptr< LayerSceneItem > LayerSceneItemHandle;

typedef std::vector< LayerSceneItemHandle > LayerGroupSceneItem;
typedef boost::shared_ptr< LayerGroupSceneItem > LayerGroupSceneItemHandle;

typedef std::vector< LayerGroupSceneItemHandle > LayerScene;
typedef boost::shared_ptr< LayerScene > LayerSceneHandle;

class LayerSceneItem
{
public:
	LayerSceneItem() {}
	~LayerSceneItem() {}

	virtual Utils::VolumeType type() = 0;

public:
	std::string layer_id_;
	double opacity_;
	Utils::GridTransform grid_transform_;
};

class DataLayerSceneItem : public LayerSceneItem
{
public:
	DataLayerSceneItem() {}
	~DataLayerSceneItem() {}

	virtual Utils::VolumeType type()
	{
		return Utils::VolumeType::DATA_E;
	}

public:
	Utils::DataVolumeHandle data_volume_;
	double contrast_;
	double brightness_;
	bool volume_rendered_;
};

class MaskLayerSceneItem : public LayerSceneItem
{
public:
	MaskLayerSceneItem() {}
	~MaskLayerSceneItem() {}

	virtual Utils::VolumeType type()
	{
		return Utils::VolumeType::MASK_E;
	}

public:
	Utils::MaskVolumeHandle mask_volume_;
	int color_;
	std::string border_;
	std::string fill_;
	bool show_isosurface_;
};

} // end namespace Seg3D

#endif