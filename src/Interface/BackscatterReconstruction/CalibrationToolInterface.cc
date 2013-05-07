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
#include <Interface/BackscatterReconstruction/CalibrationToolInterface.h>
#include <Interface/Application/LayerImporterWidget.h>

#include "ui_CalibrationToolInterface.h"

// Application Includes
#include <Application/Layer/LayerManager.h>
#include <Application/Layer/LayerGroup.h>

#include <Application/LayerIO/LayerIO.h>
//#include <Application/ViewerManager/Actions/ActionPickPoint.h>
#include <Application/PreferencesManager/PreferencesManager.h>
#include <Application/ProjectManager/ProjectManager.h>

#include <Application/Layer/Actions/ActionActivateLayer.h>

#include <Application/BackscatterReconstruction/CalibrationTool.h>


// Core Includes
#include <Core/State/Actions/ActionSetAt.h>

// test
#include <iostream>
// test

SCI_REGISTER_TOOLINTERFACE( Seg3D, CalibrationToolInterface )

namespace bfs=boost::filesystem;

namespace Seg3D
{

class OpenHISFileToolInterfacePrivate
{
public:
  //Ui::CalibrationToolInterface ui_;
  CalibrationToolInterface* interface_;
  
  void importImageStack();
  void importDataNrrd();
  void importLabelNrrd();
};

void OpenHISFileToolInterfacePrivate::importImageStack()
{
	bfs::path current_file_folder = ProjectManager::Instance()->get_current_file_folder();
  
	// Bring up the file dialog
	QString filtername;
	QStringList file_list;
  
  // Get the importer list from the LayerIO system
	std::vector< std::string > importer_types = LayerIO::Instance()->get_file_series_importer_types();
  
	QString filters;
	for ( size_t j = 0; j < importer_types.size(); j++ )
	{
    // TODO: .his and .tif only
		if ( j == 0 )
		{
			filters = QString::fromStdString( importer_types[ j ] );
			continue;
		}
		filters = filters + ";;" + QString::fromStdString( importer_types[ j ] );
	}
  
  // Use native dialog
  // TODO: not likely to be needed on Linux, but if so, test this!!!
  file_list = QFileDialog::getOpenFileNames( 0, "Select image from stack... ",
                                            current_file_folder.string().c_str(), filters, &filtername );
  
  // TODO: basically copied from LayerIOFunctions...
  // Make common utility function or just use .his and ITK .tif importers...
	if( file_list.size() == 0 )
	{
    // log no files selected
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
		QMessageBox message_box;
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
    
		QMessageBox message_box;
		message_box.setWindowTitle( "Import Layer Error" );
		message_box.addButton( QMessageBox::Ok );
		message_box.setIcon( QMessageBox::Critical );
		message_box.setText( QString::fromStdString( error_message ) );
		message_box.exec();
    QString fileStem( bfs::path( file_list.at( 0 ).toStdString() ).stem().string().c_str() );
		return;
	}
	else
	{
		// If we are able to create an importer, we then need to figure out which files are part
		// of the series.
		
		// We add the importer to the importers vector, in this case we only have a single
		// importer in the vector.
		std::vector< LayerImporterHandle > importers( 1, importer );
    
		LayerImporterWidget layer_import_dialog( importers, 0 );
		layer_import_dialog.exec();
	}  
}
  
void OpenHISFileToolInterfacePrivate::importLabelNrrd()
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

// for testing
void OpenHISFileToolInterfacePrivate::importDataNrrd()
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
  
class CalibrationToolInterfacePrivate
{
public:
	Ui::CalibrationToolInterface ui_;
	CalibrationToolInterface* interface_;  
};

// constructor
CalibrationToolInterface::CalibrationToolInterface() :
private_( new CalibrationToolInterfacePrivate ),
file_private_( new OpenHISFileToolInterfacePrivate )
{
  this->private_->interface_ = this;
  this->file_private_->interface_ = this;
}

// destructor
CalibrationToolInterface::~CalibrationToolInterface()
{
  this->disconnect_all();
}

void CalibrationToolInterface::triggerDataImport()
{
  //this->file_private_->importImageStack();
  this->file_private_->importDataNrrd();
}
  
void CalibrationToolInterface::triggerLabelImport()
{
  this->file_private_->importLabelNrrd();
}

void CalibrationToolInterface::trigger_itemActivated( QTableWidgetItem *item )
{
  LayerHandle activeLayer = LayerManager::Instance()->get_active_layer();
  if (! activeLayer)
  {
    CORE_LOG_DEBUG("No active layer");
    return;
  }

  LayerHandle layerHandle = LayerManager::Instance()->find_layer_by_name(activeLayer->get_layer_name(), -1);
  if (! layerHandle)
  {
    CORE_LOG_MESSAGE("Could not get layer by name.");
    return;
  }

  //LayerGroupHandle layerGroupHandle = layerHandle->get_layer_group();
  //std::cerr << "layerid=" << layerHandle->get_layer_id() << ", " << layerGroupHandle->
  
  QTableWidget *table = item->tableWidget();
  QList<QTableWidgetItem*> items = table->findItems( tr(activeLayer->get_layer_name().c_str()),
                                                    Qt::MatchFixedString|Qt::MatchCaseSensitive );
  if (items.size() > 0)
  {
    for (int i = 0; i < items.size(); ++i)
    {
      QTableWidgetItem* listItem = items[i];
      int listItemRow = listItem->row();
      // reset other items
      QTableWidgetItem *layerItem = table->item(listItemRow, 1);
      layerItem->setText( tr("") );
      QTableWidgetItem *colorItem = table->item(listItemRow, 2);
      colorItem->setBackground(Qt::white);
    }
  }

  int row = item->row();
  QTableWidgetItem *layerItem = table->item(row, 1);
  layerItem->setText( tr(activeLayer->get_layer_name().c_str()) );

  int color_index =  dynamic_cast< MaskLayer* >( layerHandle.get() )->color_state_->get();
  Core::Color color = PreferencesManager::Instance()->color_states_[ color_index ]->get();
  QTableWidgetItem *colorItem = table->item(row, 2);
  colorItem->setBackground(QBrush(QColor::fromRgb(color.r(), color.g(), color.b())));

  // Update tool state
	ToolHandle base_tool_ = tool();
	CalibrationTool* tool = dynamic_cast< CalibrationTool* > ( base_tool_.get() );
  tool->save( Core::Interface::GetWidgetActionContext(), row, activeLayer->get_layer_id() );
}

bool CalibrationToolInterface::build_widget( QFrame* frame )
{
  this->private_->ui_.setupUi( frame );
//	this->private_->ui_.horizontalLayout_5->setAlignment( Qt::AlignHCenter );

  this->private_->ui_.labelTableWidget->setRowCount(4);
  this->private_->ui_.labelTableWidget->setColumnCount(3);
  QStringList tableHeader;
  tableHeader << "ID" << "Label Name" << "Mask Layer";
  this->private_->ui_.labelTableWidget->setHorizontalHeaderLabels(tableHeader);
  this->private_->ui_.labelTableWidget->verticalHeader()->setVisible(false);
  this->private_->ui_.labelTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
  this->private_->ui_.labelTableWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
  this->private_->ui_.labelTableWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  this->private_->ui_.labelTableWidget->setSelectionMode(QAbstractItemView::NoSelection);
  this->private_->ui_.labelTableWidget->setShowGrid(false);

  QTableWidgetItem *index1 = new QTableWidgetItem(tr("1"));
  this->private_->ui_.labelTableWidget->setItem(0, 0, index1);
  QTableWidgetItem *layer1 = new QTableWidgetItem();
  this->private_->ui_.labelTableWidget->setItem(0, 1, layer1);
  QTableWidgetItem *color1 = new QTableWidgetItem();
  this->private_->ui_.labelTableWidget->setItem(0, 2, color1);
  QTableWidgetItem *index2 = new QTableWidgetItem(tr("2"));
  this->private_->ui_.labelTableWidget->setItem(1, 0, index2);
  QTableWidgetItem *layer2 = new QTableWidgetItem();
  this->private_->ui_.labelTableWidget->setItem(1, 1, layer2);
  QTableWidgetItem *color2 = new QTableWidgetItem();
  this->private_->ui_.labelTableWidget->setItem(1, 2, color2);
  QTableWidgetItem *index3 = new QTableWidgetItem(tr("3"));
  this->private_->ui_.labelTableWidget->setItem(2, 0, index3);
  QTableWidgetItem *layer3 = new QTableWidgetItem();
  this->private_->ui_.labelTableWidget->setItem(2, 1, layer3);
  QTableWidgetItem *color3 = new QTableWidgetItem();
  this->private_->ui_.labelTableWidget->setItem(2, 2, color3);
//  QTableWidgetItem *index4 = new QTableWidgetItem(tr("4"));
//  this->private_->ui_.labelTableWidget->setItem(3, 0, index4);
//  QTableWidgetItem *layer4 = new QTableWidgetItem();
//  this->private_->ui_.labelTableWidget->setItem(3, 1, layer4);
//  QTableWidgetItem *color4 = new QTableWidgetItem();
//  this->private_->ui_.labelTableWidget->setItem(3, 2, color4);
  
	//Step 2 - get a pointer to the tool
	ToolHandle base_tool_ = tool();
	CalibrationTool* tool = dynamic_cast< CalibrationTool* > ( base_tool_.get() );
  
	//Step 3 - connect the gui to the tool through the QtBridge
//	QtUtils::QtBridge::Connect( this->private_->ui_.input_a_, tool->target_layer_state_ );
//	QtUtils::QtBridge::Connect( this->private_->ui_.use_active_layer_, tool->use_active_layer_state_ );
//
//	QtUtils::QtBridge::Connect( this->private_->ui_.input_b_, tool->input_b_state_ );
//	QtUtils::QtBridge::Connect( this->private_->ui_.input_c_, tool->input_c_state_ );
//	QtUtils::QtBridge::Connect( this->private_->ui_.input_d_, tool->input_d_state_ );

//	QtUtils::QtBridge::Show( this->private_->ui_.message_alert_, tool->valid_target_state_, true );

//	{
//		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );	
//		this->private_->ui_.input_a_->setDisabled( tool->use_active_layer_state_->get() );
//    
//		this->connect( this->private_->ui_.use_active_layer_, SIGNAL( toggled( bool ) ),
//                  this->private_->ui_.input_a_, SLOT( setDisabled( bool ) ) );
//	}

  //connect( this->private_->ui_.input_a_, SIGNAL( currentIndexChanged(const QString & text) ), this, SLOT( trigger_table_update_a(const QString &) ) );
	QtUtils::QtBridge::Connect( this->private_->ui_.runFilterButton, boost::bind(
    &Tool::execute, tool, Core::Interface::GetWidgetActionContext() ) );
  
  connect(this->private_->ui_.labelTableWidget, SIGNAL( itemClicked( QTableWidgetItem * ) ),
    this, SLOT( trigger_itemActivated( QTableWidgetItem * ) ) );

  // TODO: make tool and action for this for scripting (take filepath, filter)
  connect( this->private_->ui_.openDataButton, SIGNAL( clicked() ), this, SLOT( triggerDataImport() ) );
  
  // TODO: make tool and action for this for scripting (take filepath, filter)
  connect( this->private_->ui_.openLabelsButton, SIGNAL( clicked() ), this, SLOT( triggerLabelImport() ) );
  
  return true;
}

} // end namespace Seg3D
