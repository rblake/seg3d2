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

// boost includes
#include<boost/tokenizer.hpp>

// Core Includes
#include <Core/Utils/Log.h>
#include <Core/Utils/StringUtil.h>
#include <Core/Interface/Interface.h>

//  Application includes
#include <Application/Tool/ToolFactory.h>
#include <Application/ToolManager/ToolManager.h>

// QtUtils includes
#include <QtUtils/Utils/QtPointer.h>
#include <QtUtils/Bridge/QtBridge.h>

// Interface Includes
#include <Interface/AppInterface/AppShortcuts.h>

#include "ui_AppShortcuts.h"

namespace Seg3D
{

class AppShortcutsPrivate
{
public:

	Ui::AppShortcuts ui_;
	QWidget* new_active_tool_;

};

AppShortcuts::AppShortcuts( QWidget *parent ) :
	QDialog( parent ),
	private_( new AppShortcutsPrivate )	
{
	// Set up the private internals of the MessageWindow class
	this->private_->ui_.setupUi( this );
	
	QIcon icon = windowIcon();
	Qt::WindowFlags flags = windowFlags();
	Qt::WindowFlags helpFlag = Qt::WindowContextHelpButtonHint;
	flags = flags & ( ~helpFlag );
	this->setWindowFlags( flags );
	this->setWindowIcon( icon );
	
	this->private_->ui_.verticalLayout_5->setAlignment( Qt::AlignTop );
	this->private_->ui_.verticalLayout_6->setAlignment( Qt::AlignTop );
	this->private_->ui_.verticalLayout_13->setAlignment( Qt::AlignTop );
	
	this->private_->ui_.active_tool_shortcuts_widget_->deleteLater();
	
	this->private_->new_active_tool_ = new QWidget();
	this->add_tool_shortcuts();

	this->add_connection( ToolManager::Instance()->activate_tool_signal_.connect( 
		boost::bind( &AppShortcuts::ShowActiveToolControls, 
		QPointer< AppShortcuts >( this ), _1 ) ) );
}

AppShortcuts::~AppShortcuts()
{
	this->disconnect_all();
}

void AppShortcuts::add_tool_shortcuts()
{
	ToolMenuList menu_list;

	ToolFactory::Instance()->list_menus( menu_list );
	ToolMenuList::const_iterator mit = menu_list.begin();
	ToolMenuList::const_iterator mit_end = menu_list.end();
	
	while ( mit != mit_end )
	{
		// Add the header widget with the name of the Tool
		QWidget* tool_shortcut_widget_ = new QWidget( this->private_->ui_.scrollAreaWidgetContents_2 );
		QSizePolicy sizePolicy1( QSizePolicy::Preferred, QSizePolicy::Fixed );
		sizePolicy1.setHorizontalStretch( 0 );
		sizePolicy1.setVerticalStretch( 0 );
		sizePolicy1.setHeightForWidth( tool_shortcut_widget_->sizePolicy(  ).hasHeightForWidth(  ) );
		tool_shortcut_widget_->setSizePolicy( sizePolicy1 );
		QVBoxLayout* verticalLayout_7 = new QVBoxLayout( tool_shortcut_widget_ );
		verticalLayout_7->setSpacing( 0 );
		verticalLayout_7->setContentsMargins( 4, 2, 4, 4 );
		verticalLayout_7->setAlignment( Qt::AlignTop );
		QWidget* tool_name_widget_ = new QWidget( tool_shortcut_widget_ );
		QHBoxLayout* horizontalLayout_33 = new QHBoxLayout( tool_name_widget_ );
		horizontalLayout_33->setSpacing( 0 );
		horizontalLayout_33->setContentsMargins( 0, 0, 0, 0 );
		QLabel* tool_name_label_ = new QLabel( tool_name_widget_ );
		tool_name_label_->setMinimumSize( QSize( 0, 20 ) );
		tool_name_label_->setText( QString::fromStdString( *mit ) + QString::fromUtf8( " -" ) );
		horizontalLayout_33->addWidget( tool_name_label_ );
		verticalLayout_7->addWidget( tool_name_widget_ );
			
		ToolInfoList tool_types_list;

		ToolFactory::Instance()->list_tools( tool_types_list, *mit );
		ToolInfoList::const_iterator it = tool_types_list.begin();
		ToolInfoList::const_iterator it_end = tool_types_list.end();
		
		while ( it != it_end )
		{
			QWidget* tool_hotkey_widget_ = new QWidget( tool_shortcut_widget_ );
			tool_hotkey_widget_->setMinimumSize( QSize( 0, 20 ) );
			QHBoxLayout* horizontalLayout_34 = new QHBoxLayout( tool_hotkey_widget_ );
			horizontalLayout_34->setSpacing( 4 );
			horizontalLayout_34->setContentsMargins( 10, 0, 0, 0 );
			QLabel* hotkey_action_ = new QLabel( tool_hotkey_widget_ );
			hotkey_action_->setText(QString::fromStdString( ( *it )->get_menu_label() ) );
			horizontalLayout_34->addWidget( hotkey_action_ );
			QLabel* hotkey_combo_ = new QLabel( tool_hotkey_widget_ );
			hotkey_combo_->setText( QString::fromStdString( ( *it )->get_shortcut_key() ) );
			horizontalLayout_34->addWidget( hotkey_combo_ );
			horizontalLayout_34->setStretch( 0, 2 );
			horizontalLayout_34->setStretch( 1, 1 );
			verticalLayout_7->addWidget( tool_hotkey_widget_ );
			++it;
		}
		
		// Add the new widget to the UI
		this->private_->ui_.verticalLayout_6->addWidget( tool_shortcut_widget_ );

		++mit;

		if( mit != mit_end )
		{
			// Add a dividing line
			QFrame* divider_line = new QFrame( this->private_->ui_.scrollAreaWidgetContents_2 );
			divider_line->setObjectName( QString::fromUtf8( "line_4" ) );
			divider_line->setFrameShape( QFrame::HLine );
			divider_line->setFrameShadow( QFrame::Sunken );
			this->private_->ui_.verticalLayout_6->addWidget( divider_line );
		}

	}

}

void AppShortcuts::ShowActiveToolControls( QPointer< AppShortcuts > qpointer, ToolHandle tool )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind(
		&AppShortcuts::show_active_tool_contols, qpointer.data(), tool ) ) );
}

void AppShortcuts::show_active_tool_contols( ToolHandle tool )
{
	std::string hotkeys = tool->get_hotkeys_and_descriptions();
	
	this->private_->new_active_tool_->deleteLater();
	
	this->private_->new_active_tool_ = new QWidget();
	this->private_->new_active_tool_->setObjectName( 
		QString::fromUtf8( "new_active_tool" ) );
	this->private_->new_active_tool_->setGeometry(QRect( 0, 0, 377, 152 ) );
	QVBoxLayout* new_verticalLayout_8 = new QVBoxLayout( this->private_->new_active_tool_ );
	new_verticalLayout_8->setSpacing( 0 );
	new_verticalLayout_8->setContentsMargins( 0, 0, 0, 0 );
	new_verticalLayout_8->setObjectName( QString::fromUtf8( "new_verticalLayout_8" ) );
	new_verticalLayout_8->setAlignment( Qt::AlignTop );
	this->private_->new_active_tool_->setStyleSheet( QString::fromUtf8( "background-color: white;" ) );
	
	this->private_->ui_.scrollArea_3->setWidget( this->private_->new_active_tool_ );
	
	if( hotkeys == "" )
	{
		this->private_->ui_.active_tool_shortcuts_->setEnabled( false );
		return;
	}
	else
	{
		this->private_->ui_.active_tool_shortcuts_->setEnabled( true );
		std::string tool_name = tool->get_menu_label();
		
		QWidget* tool_shortcut_widget_ = new QWidget( this->private_->new_active_tool_ );
		QSizePolicy sizePolicy1( QSizePolicy::Preferred, QSizePolicy::Fixed );
		sizePolicy1.setHorizontalStretch( 0 );
		sizePolicy1.setVerticalStretch( 0 );
		sizePolicy1.setHeightForWidth( tool_shortcut_widget_->sizePolicy(  ).hasHeightForWidth(  ) );
		tool_shortcut_widget_->setSizePolicy( sizePolicy1 );
		QVBoxLayout* verticalLayout_7 = new QVBoxLayout( tool_shortcut_widget_ );
		verticalLayout_7->setSpacing( 0 );
		verticalLayout_7->setContentsMargins( 4, 2, 4, 4 );
		verticalLayout_7->setAlignment( Qt::AlignTop );
		QWidget* tool_name_widget_ = new QWidget( tool_shortcut_widget_ );
		QHBoxLayout* horizontalLayout_33 = new QHBoxLayout( tool_name_widget_ );
		horizontalLayout_33->setSpacing( 0 );
		horizontalLayout_33->setContentsMargins( 0, 0, 0, 0 );
		QLabel* tool_name_label_ = new QLabel( tool_name_widget_ );
		tool_name_label_->setMinimumSize( QSize( 0, 20 ) );
		tool_name_label_->setText( QString::fromStdString( tool_name ) + QString::fromUtf8( " -" ) );
		horizontalLayout_33->addWidget( tool_name_label_ );
		verticalLayout_7->addWidget( tool_name_widget_ );
		
		std::vector< std::string > all_keys_and_descriptions = Core::SplitString( hotkeys, "|" );
		
		for( size_t i = 0; i < all_keys_and_descriptions.size(); ++i )
		{
			std::vector< std::string > key_and_description = 
				Core::SplitString( all_keys_and_descriptions[ i ], "=" );
			QWidget* tool_hotkey_widget_ = new QWidget( tool_shortcut_widget_ );
			tool_hotkey_widget_->setMinimumSize( QSize( 0, 20 ) );
			QHBoxLayout* horizontalLayout_34 = new QHBoxLayout( tool_hotkey_widget_ );
			horizontalLayout_34->setSpacing( 4 );
			horizontalLayout_34->setContentsMargins( 10, 0, 0, 0 );
			QLabel* hotkey_action_ = new QLabel( tool_hotkey_widget_ );
			hotkey_action_->setText( QString::fromStdString( key_and_description[ 0 ] ) );
			horizontalLayout_34->addWidget( hotkey_action_ );
			QLabel* hotkey_combo_ = new QLabel( tool_hotkey_widget_ );
			hotkey_combo_->setText( QString::fromStdString( key_and_description[ 1 ] ) );
			horizontalLayout_34->addWidget( hotkey_combo_ );
			horizontalLayout_34->setStretch( 0, 2 );
			horizontalLayout_34->setStretch( 1, 1 );
			verticalLayout_7->addWidget( tool_hotkey_widget_ );
		}
		
		new_verticalLayout_8->addWidget( tool_shortcut_widget_ );
	}
	
		
}








} // end namespace Seg3D