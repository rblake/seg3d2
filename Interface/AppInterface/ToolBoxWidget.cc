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

#include <Interface/AppInterface/ToolBoxWidget.h>
#include <Utils/Core/Log.h>

#include <sstream>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include <Interface/QtInterface/QtBridge.h>

// Qt includes
#include <QUrl>
#include <QDesktopServices>
#include "ui_ToolBoxPageWidget.h"

namespace Seg3D
{

typedef std::vector< QWidget* > tool_page_list_type;

class ToolBoxWidgetPrivate
{
public:
	Ui::ToolBoxPageWidget ui_;
	tool_page_list_type page_list_;
};

ToolBoxWidget::ToolBoxWidget( QWidget* parent ) :
	QScrollArea( parent )
{

	{ // Prepare the icons!!
		this->active_close_icon_.addFile( QString::fromUtf8( ":/Images/CloseWhite.png" ), QSize(),
		    QIcon::Normal, QIcon::Off );
		this->inactive_close_icon_.addFile( QString::fromUtf8( ":/Images/Close.png" ), QSize(),
		    QIcon::Normal, QIcon::Off );

		this->active_help_icon_.addFile( QString::fromUtf8( ":/Images/HelpWhite.png" ), QSize(),
		    QIcon::Normal, QIcon::Off );
		this->inactive_help_icon_.addFile( QString::fromUtf8( ":/Images/Help.png" ), QSize(),
		    QIcon::Normal, QIcon::Off );
	}

	this->private_ = ToolBoxWidgetPrivateHandle( new ToolBoxWidgetPrivate );

	setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
	setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded );
	setContentsMargins( 1, 1, 1, 1 );
	setWidgetResizable( true );

	this->main_ = new QWidget( this );
	setWidget( this->main_ );

	this->main_layout_ = new QVBoxLayout( this->main_ );
	this->main_layout_->setContentsMargins( 1, 1, 1, 1 );
	this->main_layout_->setSpacing( 1 );

	this->tool_layout_ = new QVBoxLayout;
	this->main_layout_->addLayout( this->tool_layout_ );
	this->main_layout_->addStretch();

	this->main_->setLayout( this->main_layout_ );
	this->main_->setSizePolicy( QSizePolicy::MinimumExpanding, QSizePolicy::Preferred );

}

ToolBoxWidget::~ToolBoxWidget()
{

}

void ToolBoxWidget::add_tool( QWidget * tool, const QString &label,
    boost::function< void() > close_function, boost::function< void() > activate_function,
    const std::string& help_url )
{
	if ( !tool ) return;

	//create a new base QWidget
	QWidget* new_page_ = new QWidget();
	this->private_->ui_.setupUi( new_page_ );

	this->private_->ui_.url_->setText( QString::fromStdString( help_url ) );
	this->private_->ui_.url_->hide();

	this->private_->ui_.activate_button_->setText( label );

	this->private_->ui_.help_button_->setIcon( active_help_icon_ );
	this->private_->ui_.help_button_->setIconSize( QSize( 16, 16 ) );

	this->private_->ui_.close_button_->setIcon( active_close_icon_ );
	this->private_->ui_.close_button_->setIconSize( QSize( 18, 18 ) );

	// create a new widget, send new_page_ as its parent,
	//  assign its value as tool, set its name, and add it to the tool_layout_
	QWidget* tool_ = new QWidget( new_page_ );
	tool_ = tool;
	tool_->setObjectName( QString::fromUtf8( "tool_" ) );
	this->private_->ui_.tool_frame_layout_->addWidget( tool_ );

	// add the new_page_ to the tool_layout
	this->tool_layout_->addWidget( new_page_ );

	//make all the proper connections
	connect( this->private_->ui_.help_button_, SIGNAL( clicked() ), this, SLOT(
	    help_button_clicked() ) );
	QtBridge::Connect( this->private_->ui_.activate_button_, activate_function );
	QtBridge::Connect( this->private_->ui_.close_button_, close_function );

	set_active_tool( new_page_->findChild< QWidget* > ( "tool_" ) );

	this->private_->page_list_.push_back( new_page_ );

}

void ToolBoxWidget::set_active_tool( QWidget *tool )
{
	for ( size_t i = 0; i < this->private_->page_list_.size(); i++ )
	{
		// first we deactivate the inactive tools
		if ( this->private_->page_list_[ i ]->findChild< QWidget* > ( "tool_" ) != tool )
		{
			if ( !this->private_->page_list_[ i ]->findChild< QFrame* > ( "tool_frame_" )->isHidden() )
			{
				this->private_->page_list_[ i ]->findChild< QWidget* > ( "page_background_" )->setStyleSheet(
				    QString::fromUtf8(
				        " QWidget#page_background_ { background-color: rgb(220, 220, 220); }" ) );
				this->private_->page_list_[ i ]->findChild< QPushButton* > ( "activate_button_" )->setStyleSheet(
				    QString::fromUtf8( "QPushButton#activate_button_{\n"
					    "	margin-right: 7px;\n"
					    "	height: 24px;\n"
					    "	text-align: left;\n"
					    "	padding-left: 4px;\n"
					    "	color: rgb(25, 25, 25);\n"
					    "	font: normal;\n"
					    "}\n" ) );
				this->private_->page_list_[ i ]->findChild< QToolButton* > ( "close_button_" )->setIcon(
				    inactive_close_icon_ );
				this->private_->page_list_[ i ]->findChild< QToolButton* > ( "help_button_" )->setIcon(
				    inactive_help_icon_ );
				this->private_->page_list_[ i ]->findChild< QFrame* > ( "tool_frame_" )->hide();

			}
		}

		// then, we activate the active one.
		else
		{
			this->active_index_ = static_cast< int > ( i );
			this->active_tool_ = private_->page_list_[ i ]->findChild< QWidget* > ( "tool_" );

			if ( this->private_->page_list_[ i ]->findChild< QFrame* > ( "tool_frame_" )->isHidden() )
			{
				this->private_->page_list_[ i ]->findChild< QWidget* > ( "page_background_" )->setStyleSheet(
				    QString::fromUtf8(
				        "QWidget#page_background_ { background-color: rgb(255, 128, 0); }" ) );

				this->private_->page_list_[ i ]->findChild< QPushButton* > ( "activate_button_" )->setStyleSheet(
				    QString::fromUtf8( "QPushButton#activate_button_{\n"
					    "	margin-right: 7px;\n"
					    "	height: 24px;\n"
					    "	text-align: left;\n"
					    "	padding-left: 4px;\n"
					    "	color: white;\n"
					    "	font: bold;\n"
					    "}\n" ) );
				this->private_->page_list_[ i ]->findChild< QToolButton* > ( "close_button_" )->setIcon(
				    active_close_icon_ );
				this->private_->page_list_[ i ]->findChild< QToolButton* > ( "help_button_" )->setIcon(
				    active_help_icon_ );
				this->private_->page_list_[ i ]->findChild< QFrame* > ( "tool_frame_" )->show();
			}
			
		}
	}
}

int ToolBoxWidget::index_of( QWidget *tool )
{
	for ( size_t i = 0; i < this->private_->page_list_.size(); i++ )
	{
		if ( this->private_->page_list_[ i ]->findChild< QWidget* > ( "tool_" ) == tool )
		{
			return static_cast< int > ( i );
		}
	}
	return -1;
}

QWidget* ToolBoxWidget::get_tool_at( int index )
{
	return this->private_->page_list_[ index ];
} // end get_tool_at


void ToolBoxWidget::set_active_index( int index )
{
	if ( ( index < static_cast< int > ( this->private_->page_list_.size() ) ) && ( index >= 0 ) ) set_active_tool(
	    this->private_->page_list_[ index ]->findChild< QWidget* > ( "tool_" ) );

} // end set_active_index


void ToolBoxWidget::remove_tool( int index )
{
	// Find the index that corresponds to the tool
	if ( index >=  static_cast< int > (this->private_->page_list_.size()) )
	{
		return;
	}

	this->tool_layout_->removeWidget( private_->page_list_[ index ] );
	this->private_->page_list_[ index ]->deleteLater();
	this->private_->page_list_.erase( private_->page_list_.begin() + index );

	// Set the previous tool to active if the one to be deleted is active.
	if ( this->active_index_ == index )
	{
		if ( index == 0 )
		{
			set_active_index( index );
		}
		else
		{
			set_active_index( index - 1 );
		}
	}

}

void ToolBoxWidget::help_button_clicked()
{
	QToolButton *help_button = ::qobject_cast< QToolButton* >( sender() );

	for ( size_t i = 0; i < this->private_->page_list_.size(); i++ )
	{
		if ( this->private_->page_list_[ i ]->findChild< QToolButton* > ( "help_button_" )
		    == help_button )
		{
			QDesktopServices::openUrl( QUrl(
			    this->private_->page_list_[ i ]->findChild< QLabel* > ( "url_" )->text() ) );
			break;
		}
	}
}

} // end namespace Seg3D
