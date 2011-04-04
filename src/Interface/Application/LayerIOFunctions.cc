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

// STL includes
#include <algorithm>
#include <utility>
#include <vector>
#include <set>

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/algorithm/minmax_element.hpp>

// Qt includes
#include <QMessageBox>
#include <QFileDialog>
#include <QPointer>

// Core includes
#include <Core/State/Actions/ActionSet.h>

// Application includes
#include <Application/LayerIO/LayerIO.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/LayerManager/Actions/ActionExportLayer.h>
#include <Application/PreferencesManager/PreferencesManager.h>

// Interface includes
#include <Interface/Application/LayerImporterWidget.h>
#include <Interface/Application/LayerIOFunctions.h>
#include <Interface/Application/SegmentationExportWizard.h>

namespace Seg3D
{

void LayerIOFunctions::ImportFiles( QMainWindow* main_window, std::string file_to_open )
{
	// List of files that need to be imported.
	QStringList file_list;
	// The name of the importer to use.
	std::string importer_name;
	
	if( file_to_open != "" )
	{
		// If there is a file, do not open a dialog and just add it into the file_list.
		file_list << QString::fromStdString( file_to_open );
		// Let the program find the right importer.
		importer_name = "All Importers (*)";
	}
	else
	{
		// Get the importer list from the LayerIO system.
		std::vector< std::string > importer_types = 
			LayerIO::Instance()->get_single_file_importer_types();
		
		// Make a list to select the right importer from.
		QString filters = "";
		for ( size_t j = 0; j < importer_types.size(); j++ )
		{
			if( j == 0 )
			{
				filters = QString::fromStdString( importer_types[j] );
				continue;
			}
			filters = filters + ";;" + QString::fromStdString( importer_types[j] );
		}

		// Bring up the file dialog.
		QString qs_filtername;
		file_list = QFileDialog::getOpenFileNames( main_window, 
			"Import Layer(s)... ", "/home", filters, &qs_filtername );
		
		// If no file was selected just return
		if( file_list.size() == 0 )
		{
			QMessageBox message_box( main_window );
			message_box.setWindowTitle( "Import Layer Error" );
			message_box.addButton( QMessageBox::Ok );
			message_box.setIcon( QMessageBox::Critical );
			message_box.setText( "No files were selected." );
			message_box.exec();	
			return;
		}

		// Get the name of the importer that was selected
		importer_name = qs_filtername.toStdString();
	}

	// List of importers to use
	std::vector< LayerImporterHandle > importers;

	// Loop over all the selected files and generate an importer for them
	for( int i = 0; i < file_list.size(); ++i )
	{
		LayerImporterHandle importer;
		
		// Create a new importer
		if( ! ( LayerIO::Instance()->create_single_file_importer( file_list.at( i ).toStdString(), 
			importer, importer_name ) ) )
		{
			// Failed to create the importer, and warn the user explicitly
			std::string error_message = std::string( "ERROR: No importer is available for file '" ) 
				+ file_list.at( i ).toStdString() + std::string( "'." );

			QMessageBox message_box( main_window );
			message_box.setWindowTitle( "Import Layer Error" );
			message_box.addButton( QMessageBox::Ok );
			message_box.setIcon( QMessageBox::Critical );
			message_box.setText( QString::fromStdString( error_message ) );
			
			// Send message box to Qt
			message_box.exec();
			return;
		}
		else
		{
			// Add the importer to the list of all the importers that need to be handled
			importers.push_back( importer );
		}
	}

	// Open the importer dialog that issues the action to import the data file(s)
	LayerImporterWidget layer_import_dialog( importers, main_window );
	layer_import_dialog.exec();
}

void LayerIOFunctions::ImportSeries( QMainWindow* main_window )
{
	// Get the importer list from the LayerIO system
	std::vector< std::string > importer_types = LayerIO::Instance()->get_file_series_importer_types();

	QString filters = "";
	for ( size_t j = 0; j < importer_types.size(); j++ )
	{
		if( j == 0 )
		{
			filters = QString::fromStdString( importer_types[ j ] );
			continue;
		}
		filters = filters + ";;" + QString::fromStdString( importer_types[ j ] );
	}

	// Bring up the file dialog
	QString filtername;
	QStringList file_list = QFileDialog::getOpenFileNames( main_window, 
		"Select a file from the series... ", "/home", filters, &filtername );

	// If no files were selected just exit
	if( file_list.size() == 0 )
	{
		QMessageBox message_box( main_window );
		message_box.setWindowTitle( "Import Layer Error" );
		message_box.addButton( QMessageBox::Ok );
		message_box.setIcon( QMessageBox::Critical );
		message_box.setText( "No files were selected."  );
		message_box.exec();	
		return;
	}

	// Convert the filenames into STL strings
	std::vector< std::string > filenames( file_list.size() );
	for ( int j = 0; j < file_list.size(); j++ )
	{
		filenames[ j ] = file_list.at( j ).toStdString();
	}
		
	// Find all the filenames that are associated with this series
	if ( !( LayerIO::FindFileSeries( filenames ) ) )
	{
		QMessageBox message_box( main_window );
		message_box.setWindowTitle( "Import Layer Error" );
		message_box.addButton( QMessageBox::Ok );
		message_box.setIcon( QMessageBox::Critical );
		message_box.setText( "Could not resolve filenames." );
		message_box.exec();	
		return;
	}
		
	std::string importer_name = filtername.toStdString();

	LayerImporterHandle importer;
	if( ! ( LayerIO::Instance()->create_file_series_importer( filenames, importer, 
		importer_name ) ) )
	{
		// If we are unable to create an importer we pop up an error message box
		std::string error_message = std::string( "ERROR: No importer is available for file '" ) + 
			file_list.at( 0 ).toStdString() + std::string( "'." );

		QMessageBox message_box( main_window );
		message_box.setWindowTitle( "Import Layer Error" );
		message_box.addButton( QMessageBox::Ok );
		message_box.setIcon( QMessageBox::Critical );
		message_box.setText( QString::fromStdString( error_message ) );
		message_box.exec();
		return;
	}
	else
	{
		// If we are able to create an importer, we then need to figure out which files are part
		// of the series.
		
		// We add the importer to the importers vector, in this case we only have a single
		// importer in the vector.
		std::vector< LayerImporterHandle > importers( 1, importer );

		LayerImporterWidget layer_import_dialog( importers, main_window );
		layer_import_dialog.exec();
	}
}
	
void LayerIOFunctions::ExportLayer( QMainWindow* main_window )
{
	std::vector< LayerHandle > layer_handles;
	layer_handles.push_back( LayerManager::Instance()->get_active_layer() );
	if( !layer_handles[ 0 ] ) return;

	if( layer_handles[ 0 ]->get_type() != Core::VolumeType::DATA_E )
	{
		std::string error_message = 
			std::string( "ERROR: A Data layer is not set as the active layer" );

		QMessageBox message_box( main_window );
		message_box.setWindowTitle( "Export Layer Error." );
		message_box.addButton( QMessageBox::Ok );
		message_box.setIcon( QMessageBox::Critical );
		message_box.setText( QString::fromStdString( error_message ) );
		message_box.exec();
		return;
	}

	boost::filesystem::path file_path = boost::filesystem::path( 
		PreferencesManager::Instance()->export_path_state_->get() ) / layer_handles[ 0 ]->get_layer_name();

	QString filename = QFileDialog::getSaveFileName( main_window, "Export Data Layer As... ",
		QString::fromStdString( file_path.string() ),
		"NRRD files (*.nrrd);;DICOM files (*.dcm);;TIFF files (*.tiff);;PNG files (*.png)" );
	
	if( filename == "" ) return;
	
	if( boost::filesystem::exists( boost::filesystem::path( filename.toStdString() ).parent_path() ) )
	{
		Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(),
			PreferencesManager::Instance()->export_path_state_, 
			boost::filesystem::path( filename.toStdString() ).parent_path().string() );
	}
		
	std::string extension = boost::filesystem::path( filename.toStdString() ).extension(); 
	std::string exportername;
	
	if( extension == ".nrrd" ) exportername = "NRRD Exporter";
	else if( extension != "" ) exportername = "ITK Data Exporter";
		
	LayerExporterHandle exporter;
	if( ! ( LayerIO::Instance()->create_exporter( exporter, layer_handles, exportername, extension ) ) )
	{
		std::string error_message = std::string("ERROR: No exporter is available for file '") + 
			filename.toStdString() + std::string("'.");

		QMessageBox message_box( main_window );
		message_box.setWindowTitle( "Import Layer..." );
		message_box.addButton( QMessageBox::Ok );
		message_box.setIcon( QMessageBox::Critical );
		message_box.setText( QString::fromStdString( error_message ) );
		message_box.exec();
		return;
	}
		
	ActionExportLayer::Dispatch( Core::Interface::GetWidgetActionContext(), 
		layer_handles[ 0 ]->get_layer_id(), exporter, filename.toStdString() );
}

void LayerIOFunctions::ExportSegmentation( QMainWindow* main_window )
{
	QPointer< SegmentationExportWizard > export_segmentation_wizard_ = 
		new SegmentationExportWizard( main_window );
	export_segmentation_wizard_->show();
}

} // end namespace Seg3D
