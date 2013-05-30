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


// Boost includes
#include <boost/filesystem.hpp>

// Core includes
#include <Core/Interface/Interface.h>
#include <Core/Utils/Log.h>

// Qt includes
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>
#include <QtGui/QClipboard>
#include <QtGui/QColorDialog>
#include <QtGui/QBrush>
#include <QtGui/QColor>

// Qt Gui Includes
#include <QtUtils/Bridge/QtBridge.h>

// Core Includes
#include <Core/Utils/Log.h>
#include <Core/State/Actions/ActionSet.h>

// Interface Includes
#include <Interface/BackscatterReconstruction/ViewerToolInterface.h>
#include <Interface/Application/LayerImporterWidget.h>

#include "ui_ViewerToolInterface.h"

// Application Includes
#include <Application/Layer/LayerManager.h>
#include <Application/LayerIO/LayerIO.h>
#include <Application/PreferencesManager/PreferencesManager.h>
#include <Application/ProjectManager/ProjectManager.h>

#include <Application/Layer/Actions/ActionActivateLayer.h>


// Core Includes
#include <Core/State/Actions/ActionSetAt.h>

// test
#include <iostream>
// test

SCI_REGISTER_TOOLINTERFACE( Seg3D, ViewerToolInterface )

namespace bfs=boost::filesystem;

namespace Seg3D
{

class ViewerToolInterfacePrivate
{
public:
  Ui::ViewerToolInterface ui_;
  ViewerToolInterface* interface_;
  
  //void importImageStack();
  void importLabelNrrd();
};

void ViewerToolInterfacePrivate::importLabelNrrd()
{
  bfs::path current_file_folder = ProjectManager::Instance()->get_current_file_folder();
  
  // Bring up the file dialog
  QString filtername;
  QStringList file_list;
  
  // Get the importer list from the LayerIO system
  std::vector< std::string > importer_types = LayerIO::Instance()->get_single_file_importer_types();
  
  QString filters;
  for ( size_t j = 0; j < importer_types.size(); j++ )
  {
    // TODO: .nrrd only
    if ( j == 0 )
    {
      filters = QString::fromStdString( importer_types[ j ] );
      continue;
    }
    filters = filters + ";;" + QString::fromStdString( importer_types[ j ] );
  }
  
  // Use native dialog
  // TODO: not likely to be needed on Linux, but if so, test this!!!
  file_list = QFileDialog::getOpenFileNames( 0, "Select labels file ",
                                            current_file_folder.string().c_str(), filters, &filtername );
  
  // TODO: basically copied from LayerIOFunctions...
  // Make common utility function or just use .his and ITK .tif importers...
  if( file_list.size() == 0 )
  {
    // log no files selected
    return;
  }
  
  std::string importer_name = filtername.toStdString();
  
  LayerImporterHandle importer;
  // TODO: specialize to just masks?
  if( ! ( LayerIO::Instance()->create_single_file_importer( file_list.at( 0 ).toStdString(), importer, 
                                                           importer_name ) ) )
  {
    // If we are unable to create an importer we pop up an error message box
    std::string error_message = std::string( "ERROR: No importer is available for file '" ) + 
    file_list.at( 0 ).toStdString() + std::string( "'." );
    
    QMessageBox message_box;
    message_box.setWindowTitle( "Import Layer Error" );
    message_box.addButton( QMessageBox::Ok );
    message_box.setIcon( QMessageBox::Critical );
    message_box.setText( QString::fromStdString( error_message ) );
    message_box.exec();
    QString fileStem( bfs::path( file_list.at( 0 ).toStdString() ).stem().string().c_str() );
    //ui_.fileNameLabel->setText( fileStem );
    return;
  }
  else
  {
    // List of importers to use
    std::vector< LayerImporterHandle > importers;
    
    // Loop over all the selected files and generate an importer for them
    for( int i = 0; i < file_list.size(); ++i )
    {
      LayerImporterHandle importer;
      std::string error;
      
      // Create a new importer
      if( ! ( LayerIO::Instance()->create_single_file_importer( file_list.at( i ).toStdString(), 
                                                               importer, error, importer_name ) ) )
      {
        // Failed to create the importer, and warn the user explicitly
        std::string error_message = std::string( "ERROR: No importer is available for file '" ) 
        + file_list.at( i ).toStdString() + std::string( "'. " ) + error;
        
        QMessageBox message_box;
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
    LayerImporterWidget layer_import_dialog( importers, 0 );
    layer_import_dialog.exec();
  }
}

// constructor
ViewerToolInterface::ViewerToolInterface() :
private_( new ViewerToolInterfacePrivate )
{
  this->private_->interface_ = this;
}

// destructor
ViewerToolInterface::~ViewerToolInterface()
{
  this->disconnect_all();
}

//void ViewerToolInterface::triggerDataImport()
//{
//  this->private_->importImageStack();
//}

void ViewerToolInterface::triggerLabelImport()
{
  this->private_->importLabelNrrd();
}

bool ViewerToolInterface::build_widget( QFrame* frame )
{
  this->private_->ui_.setupUi( frame );

  // TODO: make tool and action for this for scripting (take filepath, filter)
  connect( this->private_->ui_.openLabelsButton, SIGNAL( clicked() ), this, SLOT( triggerLabelImport() ) );
  
  return true;
}
  
} // end namespace Seg3D
