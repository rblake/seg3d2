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

// Boost includes
#include <boost/filesystem.hpp>

// Application includes
#include <Application/UndoBuffer/UndoBuffer.h>
#include <Application/LayerIO/LayerIO.h>
#include <Application/LayerManager/Actions/ActionImportLayer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/LayerManager/LayerUndoBufferItem.h>
#include <Application/ProjectManager/ProjectManager.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
CORE_REGISTER_ACTION( Seg3D, ImportLayer )

namespace Seg3D
{

bool ActionImportLayer::validate( Core::ActionContextHandle& context )
{
	// Validate whether the filename is actually valid
	
	// Convert the file to a boost filename
	boost::filesystem::path full_filename( this->filename_ );
	
	// Check whether the file exists
	if ( !( boost::filesystem::exists( full_filename ) ) )
	{
		context->report_error( std::string( "File '" ) + this->filename_ + "' does not exist." );
		return false;
	}

	// Update the filename to include the full path, so the filename has an absolute path
	try
	{
		full_filename = boost::filesystem::absolute( full_filename );
	}
	catch ( ... )
	{
		context->report_error(  std::string( "Could not determine full path of '" ) +
			this->filename_ + "'." );
		return false;
	}
	
	// Reinsert the filename in the action
	this->filename_ = full_filename.string();
	
	// If there is no layer importer we need to generate one
	if ( !this->layer_importer_ )
	{
		// Create a new importer
		if ( !( LayerIO::Instance()->create_single_file_importer( this->filename_,  
			this->layer_importer_, this->importer_ ) ) )
		{
			// In case the importer could not be generated
			context->report_error( std::string( "Could not create importer with name '" ) +
				this->importer_ + "' for file '" + this->filename_ + "'." );
			return false;
		}
		
		if ( !this->layer_importer_ )
		{
			context->report_error( std::string( "Could not create importer with name '" ) +
				this->importer_ + "' for file '" + this->filename_ + "'." );
			return false;			
		}
	}

	// Check whether mode is a valid string
	if ( this->mode_ != "data" && this->mode_ != "single_mask" && this->mode_ != "bitplane_mask" &&
		this->mode_ != "label_mask" )
	{
		context->report_error( "Importer mode needs to be data, single_mask, bitplane_mask, or"
			" label_mask." );
		return false;
	}
	
	// Check the information that we can retrieve from the header of this file
	LayerImporterFileInfoHandle info;
	if ( !( this->layer_importer_->get_file_info( info ) ) )
	{
		context->report_error( this->layer_importer_->get_error() );
		return false;
	} 

	if ( this->mode_ != "data" && !info->get_mask_compatible() )
	{
		context->report_error( std::string( "Import mode '") +  this->mode_ + 
			"' is not available for this importer." );
		return false;
	}

	return true; // validated
}


bool ActionImportLayer::run( Core::ActionContextHandle& context, Core::ActionResultHandle& result )
{
	// Get the current counters for groups and layers, so we can undo the changes to those counters
	// NOTE: This needs to be done before a new layer is created
	LayerManager::id_count_type id_count = LayerManager::GetLayerIdCount();

	// Forwarding a message to the UI that we are importing a layer. This generates a progress bar
	std::string message = std::string( "Importing '" ) + this->filename_ + std::string( "'" );
	Core::ActionProgressHandle progress = Core::ActionProgressHandle( 
		new Core::ActionProgress( message ) );

	// Indicate that we have started the process
	progress->begin_progress_reporting();
	
	// The ImporterFileData is an abstraction of all the data can be extracted from the file
	LayerImporterFileDataHandle data;
	
	// Get the data from the file
	if ( !( this->layer_importer_->get_file_data( data ) ) )
	{
		progress->end_progress_reporting();
		context->report_error( "Layer importer failed to extract volume data from file." );
		return false;
	}
	
	// Now convert this abtract intermediate into layers that can be inserted in the program
	// NOTE: This step is only reformating the header of the data and adds the state variables
	// for the layers.
	std::vector< LayerHandle > layers;
	if ( !( data->convert_to_layers( this->mode_, layers ) ) )
	{
		progress->end_progress_reporting();
		context->report_error( "Importer could not convert data into the requested format." );
		return false;	
	}
	
	// Now insert the layers one by one into the layer manager.
	for ( size_t j = 0; j < layers.size(); j++ )
	{
		layers[ j ]->provenance_id_state_->set( this->get_output_provenance_id( j ) );
		LayerManager::Instance()->insert_layer( layers[ j ] );
	}

	// Now the layers are properly inserted, generate the undo item that will undo this action.
	{
		// Create an undo item for this action
		LayerUndoBufferItemHandle item( new LayerUndoBufferItem( "Import Layer" ) );

		// Tell which action has to be re-executed to obtain the result
		item->set_redo_action( this->shared_from_this() );

		// Tell which layer was added so undo can delete it
		for ( size_t j = 0; j < layers.size(); j++ )
		{
			item->add_layer_to_delete( layers[ j ] );
		}
		// Tell what the layer/group id counters are so we can undo those as well
		item->add_id_count_to_restore( id_count );
		
		// Add the complete record to the undo buffer
		UndoBuffer::Instance()->insert_undo_item( context, item );

		// Create a provenance record
		ProvenanceStepHandle provenance_step( new ProvenanceStep );
		
		// Get the input provenance ids from the translate step
		provenance_step->set_input_provenance_ids( this->get_input_provenance_ids() );
		
		// Get the output and replace provenance ids from the analysis above
		provenance_step->set_output_provenance_ids( this->get_output_provenance_ids() );
			
		// Get the action and turn it into provenance	
		provenance_step->set_action( this->export_to_provenance_string() );		
		
		// Add step to provenance record
		ProjectManager::Instance()->get_current_project()->add_to_provenance_database(
			provenance_step );		
	}
	
	boost::filesystem::path full_filename( this->filename_ );
	ProjectManager::Instance()->current_file_folder_state_->set( 
		full_filename.parent_path().string() );
	ProjectManager::Instance()->checkpoint_projectmanager();

	// We are done processing
	progress->end_progress_reporting();

	return true;
}

void ActionImportLayer::clear_cache()
{
	// After the action has been executed this one needs to be reset, as it should not linger around
	// because it has the data volume associated with it. Hence to conserve memory we will 
	// regenerate this one when the action is executed again.
	this->layer_importer_.reset();
}

void ActionImportLayer::Dispatch( Core::ActionContextHandle context, const std::string& filename, 
	const std::string& mode, const std::string importer )
{
	// Create new action
	ActionImportLayer* action = new ActionImportLayer;
	
	// Set action parameters
	action->filename_ = filename;
	action->importer_ = importer;
	action->mode_ = mode;
	action->cache_ = -1;
	
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

void ActionImportLayer::Dispatch( Core::ActionContextHandle context, 
	const LayerImporterHandle& importer, const std::string& mode )
{
	// Create new action
	ActionImportLayer* action = new ActionImportLayer;
	
	// Fill in the short cut so the data that was already read is not lost
	action->layer_importer_ = importer;

	// We need to fill in these to ensure the action can be replayed without the importer present
	action->filename_ = importer->get_filename();
	action->importer_ = importer->get_name();
	action->mode_ = mode;
	action->cache_ = -1;
	
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}
	
} // end namespace Seg3D
