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

#ifndef APPLICATION_LAYERMANAGER_ACTIONS_ACTIONIMPORTSERIES_H
#define APPLICATION_LAYERMANAGER_ACTIONS_ACTIONIMPORTSERIES_H

// Core includes
#include <Core/Action/Actions.h>
#include <Core/Interface/Interface.h>

// Application includes
#include <Application/LayerIO/LayerImporter.h>

namespace Seg3D
{

// TODO: We should split this in importfromfile and importfromseries
// --JGS
	
class ActionImportSeries : public Core::Action
{

CORE_ACTION( 
	CORE_ACTION_TYPE( "ImportSeries", "This action imports a series into the layer manager.")
	CORE_ACTION_ARGUMENT( "filename", "The name of the file to load." )
	CORE_ACTION_KEY( "mode", "data", "The mode to use: data, single_mask, bitplane_mask, or label_mask.")
	CORE_ACTION_KEY( "importer", "", "Optional name for a specific importer." )
	
	CORE_ACTION_CHANGES_PROJECT_DATA()
)

	// -- Constructor/Destructor --
public:
	ActionImportSeries()
	{
		this->add_argument( this->filename_ );
		this->add_key( this->mode_ );
		this->add_key( this->importer_ );
	}
	
	virtual ~ActionImportSeries()
	{
	}
	
	// -- Functions that describe action --
public:
	// VALIDATE:
	// Each action needs to be validated just before it is posted. This way we
	// enforce that every action that hits the main post_action signal will be
	// a valid action to execute.
	virtual bool validate( Core::ActionContextHandle& context );

	// RUN:
	// Each action needs to have this piece implemented. It spells out how the
	// action is run. It returns whether the action was successful or not.
	virtual bool run( Core::ActionContextHandle& context, Core::ActionResultHandle& result );

	// CLEAR_CACHE:
	// Clear any objects that were given as a short cut to improve performance.
	virtual void clear_cache();	
		
	// -- Action parameters --
private:

	// The filename of the file to load
	Core::ActionParameter< std::string > filename_;

	// How should the file be loaded
	Core::ActionParameter< std::string > mode_;

	// Which type of importer should we use
	Core::ActionParameter< std::string > importer_;

	// Short cut to the layer importer that has already loaded the data if the file
	// was read through the GUI
	LayerImporterHandle layer_importer_;
	
	// -- Dispatch this action from the interface --
public:
	// CREATE:
	// Create action that imports a layer
	static Core::ActionHandle Create( const std::string& filename, const std::string& mode = "data",
		const std::string importer = "" );

	// CREATE:
	// Create action that imports a layer
	static Core::ActionHandle Create( const LayerImporterHandle& importer, LayerImporterMode mode );
	
	// DISPATCH:
	// Create and dispatch action that moves the layer above 
	static void Dispatch( Core::ActionContextHandle context, const std::string& filename, 
		const std::string& mode = "data", const std::string importer = "" );

	// DISPATCH:
	// To avoid reading a file twice, this action has a special option, so it can take an
	// importer that has already loaded the file. This prevents it from being read twice
	static void Dispatch( Core::ActionContextHandle context, const LayerImporterHandle& importer, 
		LayerImporterMode mode );
	
};
	
} // end namespace Seg3D

#endif
