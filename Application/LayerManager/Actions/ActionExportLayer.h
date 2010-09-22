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

#ifndef APPLICATION_LAYERMANAGER_ACTIONS_ACTIONEXPORTLAYER_H
#define APPLICATION_LAYERMANAGER_ACTIONS_ACTIONEXPORTLAYER_H

// Core includes
#include <Core/Action/Actions.h>
#include <Core/Interface/Interface.h>

#include <Application/LayerIO/LayerExporter.h>

namespace Seg3D
{
	
class ActionExportLayer : public Core::Action
{

CORE_ACTION( 
	CORE_ACTION_TYPE( "ExportLayer", "This action exports a layer to file.")
	CORE_ACTION_ARGUMENT( "file_path", "A path, including the name of the file where the layer should be exported to." )
	CORE_ACTION_KEY( "mode", "data", "The mode to use: data, single_mask, bitplane_mask, or label_mask.")
	CORE_ACTION_KEY( "exporter", "", "Optional name for a specific exporter." )
)

	// -- Constructor/Destructor --
public:
	ActionExportLayer()
	{
		add_argument( this->file_path_ );
		add_key( this->mode_ );
		add_key( this->exporter_ );
	}
	
	virtual ~ActionExportLayer()
	{
	}
	
	// -- Functions that describe action --
public:
	virtual bool validate( Core::ActionContextHandle& context );
	virtual bool run( Core::ActionContextHandle& context, Core::ActionResultHandle& result );
	
	// -- Action parameters --
private:
	// Where the layer should be exported
	Core::ActionParameter< std::string > file_path_;
	
	// How should the file be exported
	Core::ActionParameter< std::string > mode_;
	
	// Which type of exporter should we use
	Core::ActionParameter< std::string > exporter_;
	
	// Short cut to the layer exporter that has already loaded the data if the file
	// was read through the GUI
	Core::ActionCachedHandle<LayerExporterHandle> layer_exporter_;
	
	// -- Dispatch this action from the interface --
public:
	// CREATE:
	// Create action that exports a segmentation
	static Core::ActionHandle Create( const LayerExporterHandle& exporter, LayerExporterMode mode,
		const std::string& file_path );

	// DISPATCH:
	// To avoid reading a file twice, this action has a special option, so it can take an
	// importer that has already loaded the file. This prevents it from being read twice
	static void Dispatch( Core::ActionContextHandle context, const LayerExporterHandle& exporter, 
		LayerExporterMode mode, const std::string& file_path );
	
};
	
} // end namespace Seg3D

#endif
