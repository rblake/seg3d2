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

// STD includes
#include <sstream>
#include <iostream>

// Boost includes
#include <boost/lexical_cast.hpp>

// Qt Includes
#include <QVBoxLayout>

// Core includes
#include <Core/Utils/Log.h>

// QtUtils includes
#include <QtUtils/Utils/QtPointer.h>

//Application Includes
#include <Application/Layer/LayerManager.h>

//Interface Includes
#include <Interface/Application/LayerManagerWidget.h>
#include <Interface/Application/LayerGroupWidget.h>

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// Class LayerManagerWidgetPrivate
//////////////////////////////////////////////////////////////////////////

typedef std::map< std::string, LayerGroupWidget* > GroupWidgetMap;

class LayerManagerWidgetPrivate : public QObject
{
	// Constructor
public:
	LayerManagerWidgetPrivate( LayerManagerWidget* parent ) :
		QObject( parent ),
		parent_( parent ),
		loading_states_( false )
	{
	}

	// -- functions that handle signals --
public:
	// HANDLE_LAYER_INSERTED:
	// This function is called in response to the LayerManager::layer_inserted_signal_.
	void handle_layer_inserted( LayerHandle layer, bool new_group );

	// HANDLE_LAYERS_DELETED:
	// This function is called in response to the LayerManager::layers_deleted_signal_.
	void handle_layers_deleted( std::vector< std::string > affected_groups, bool groups_deleted );

	// HANDLE_LAYERS_REORDERED:
	// This function is called in response to LayerManager::layers_reordered_signal_.
	void handle_layers_reordered( std::string group_id );

	// HANDLE_GROUPS_REORDERED:
	// This function is called in response to LayerManager::groups_reordered_signal_.
	void handle_groups_reordered();

	// UPDATE_GROUP_WIDGETS:
	// Update the UI to reflect the current data in LayerManager.
	void update_group_widgets();

	// PRE_LOAD_STATES:
	// Called on pre_load_states_signal_ from the StateEngine.
	void pre_load_states();

	// POST_LOAD_STATES:
	// Called on post_load_states_signal_ from the StateEngine.
	void post_load_states();

	// RESET:
	// Remove all the layer widgets. Called on Core::Application::reset_signal_.
	void reset();

	// MAKE_NEW_GROUP:
	// function that creates a new group to put layers into. 
	LayerGroupWidget* make_new_group( LayerGroupHandle group );

public:
	// Pointer to the LayerManagerWidget
	LayerManagerWidget* parent_;
	// Main widget inside the scroll area
	QWidget* main_;
	
	// Layouts that define where the groups are located within the widget
	QVBoxLayout* main_layout_;
	QVBoxLayout* group_layout_;

	// Mapping of group name to its underlying widget
	GroupWidgetMap group_map_;
	// A flag indicating a session loading in process
	bool loading_states_;

	// -- static functions for callbacks into this widget --
public:
	typedef QPointer< LayerManagerWidgetPrivate > qpointer_type;

	// HANDLELAYERINSERTED:
	// This function is called in response to the LayerManager::layer_inserted_signal_.
	static void HandleLayerInserted( qpointer_type qpointer, LayerHandle layer, bool new_group );

	// HANDLELAYERSDELETED:
	// This function is called in response to the LayerManager::layers_deleted_signal_.
	static void HandleLayersDeleted( qpointer_type qpointer, 
		std::vector< std::string > affected_groups, bool groups_deleted );

	// HANDLELAYERSREORDERED:
	// This function is called in response to LayerManager::layers_reordered_signal_.
	static void HandleLayersReordered( qpointer_type qpointer, std::string group_id );

	// HANDLEGROUPSREORDERED:
	// This function is called in response to LayerManager::groups_reordered_signal_.
	static void HandleGroupsReordered( qpointer_type qpointer );

	// PRELOADSTATES:
	// Called on pre_load_states_signal_ from the StateEngine.
	static void PreLoadStates( qpointer_type qpointer );

	// POSTLOADSTATES:
	// Called on post_load_states_signal_ from the StateEngine.
	static void PostLoadStates( qpointer_type qpointer );

	// RESET:
	// Called on Core::Application::reset_signal_.
	static void HandleReset( qpointer_type qpointer );
};

void LayerManagerWidgetPrivate::handle_layer_inserted( LayerHandle layer, bool new_group )
{
	// Don't update if there is a session loading going on
	if ( this->loading_states_ ) return;

	LayerGroupHandle layer_group = layer->get_layer_group();
	// The layer group no longer exists, ignore.
	// NOTE: This could happen if a layer is created and then deleted very quickly.
	// It's not likely to happen through the GUI, but very likely through scripting.
	if ( !layer_group )
	{
		return;
	}

	// Disable UI updates
	this->parent_->setUpdatesEnabled( false );

	// If the inserted layer belongs to a new group, we need to update the group widgets
	if ( new_group )
	{
		this->update_group_widgets();
	}
	else
	{
		GroupWidgetMap::iterator it = this->group_map_.find( layer_group->get_group_id() );
		// Make sure that the layer group still exists
		if ( it != this->group_map_.end() )
		{
			it->second->handle_change();
		}
	}

	// Enable UI updates to redraw
	this->parent_->setUpdatesEnabled( true );
}

void LayerManagerWidgetPrivate::handle_layers_deleted( 
	std::vector< std::string > affected_groups, bool groups_deleted )
{
	// Don't update if there is a session loading going on
	if ( this->loading_states_ ) return;

	// Disable UI updates
	this->parent_->setUpdatesEnabled( false );
	
	// Update group widgets if some groups have been deleted
	if ( groups_deleted )
	{
		this->update_group_widgets();
	}

	// For each affected group, update the group widget
	for ( size_t i = 0; i < affected_groups.size(); ++i )
	{
		GroupWidgetMap::iterator widget_it = this->group_map_.find( affected_groups[ i ] );
		if ( widget_it != this->group_map_.end() )
		{
			widget_it->second->handle_change();
		}
	}

	// Enable UI updates to redraw
	this->parent_->setUpdatesEnabled( true );
}

void LayerManagerWidgetPrivate::handle_layers_reordered( std::string group_id )
{
	// Don't update if there is a session loading going on
	if ( this->loading_states_ ) return;

	GroupWidgetMap::iterator widget_it = this->group_map_.find( group_id );
	if ( widget_it != this->group_map_.end() )
	{
		widget_it->second->handle_change();
	}
}

void LayerManagerWidgetPrivate::handle_groups_reordered()
{
	// Don't update if there is a session loading going on
	if ( this->loading_states_ ) return;

	this->parent_->setUpdatesEnabled( false );
	this->update_group_widgets();
	this->parent_->setUpdatesEnabled( true );
}

void LayerManagerWidgetPrivate::update_group_widgets()
{
	// Don't update if there is a session loading going on
	if ( this->loading_states_ ) return;

	// Get a list of all the groups
	std::vector< LayerGroupHandle > groups;
	LayerManager::Instance()->get_groups( groups );
	// Make a copy of the old widgets map
	GroupWidgetMap tmp_map = this->group_map_;
	// Clear the original map
	this->group_map_.clear();

	// Loop through all the current groups and put their corresponding widgets
	// in the right order. Create new widgets when necessary.
	for ( size_t i = 0; i < groups.size(); ++i )
	{
		// Look for an existing widget for the group. If found, remove it from the 
		// temporary map. Otherwise, create a new widget.
		LayerGroupWidget* group_widget;
		std::string group_id = groups[ i ]->get_group_id();
		GroupWidgetMap::iterator it = tmp_map.find( group_id );
		if ( it != tmp_map.end() )
		{
			group_widget = it->second;
			tmp_map.erase( it );
			
			// Instantly hide any still visible drop space
			group_widget->instant_hide_drop_space();
		}
		else
		{
			group_widget = this->make_new_group( groups[ i ] );
			group_widget->handle_change();
			group_widget->show();
		}

		// Put the widget in the layout
		this->group_layout_->insertWidget( static_cast< int >( i ), group_widget );
		// Add the widget to the map
		this->group_map_[ group_id ] = group_widget;
	}

	// For anything left in the temporary map, they are no longer needed.
	// Remove them from the layout and delete them.
	for ( GroupWidgetMap::iterator it = tmp_map.begin(); it != tmp_map.end(); ++it )
	{
		LayerGroupWidget* group_widget = it->second;
		// Remove the widget from the layout and mark it for deletion
		this->group_layout_->removeWidget( group_widget );
		group_widget->deleteLater();
	}
	tmp_map.clear();

	// Squeeze the LayerManagerWidget
	this->parent_->setMinimumHeight( 0 );
}

void LayerManagerWidgetPrivate::pre_load_states()
{
	// All the widgets should have been deleted before session loading started
	assert( this->group_map_.empty() );

	this->loading_states_ = true;
}

void LayerManagerWidgetPrivate::post_load_states()
{
	this->loading_states_ = false;
	this->parent_->setUpdatesEnabled( false );
	this->update_group_widgets();
	this->parent_->setUpdatesEnabled( true );
}

LayerGroupWidget* LayerManagerWidgetPrivate::make_new_group( LayerGroupHandle group )
{
	LayerGroupWidget* new_group = new LayerGroupWidget( this->parent_, group );

	connect( new_group, SIGNAL( prep_layers_for_drag_and_drop_signal_( bool ) ), 
		this->parent_, SLOT( prep_layers_for_drag_and_drop( bool ) ) );

	connect( new_group, SIGNAL( prep_groups_for_drag_and_drop_signal_( bool ) ), 
		this->parent_, SLOT( prep_groups_for_drag_and_drop( bool ) ) );

	connect( new_group, SIGNAL( picked_up_group_size_signal_( int ) ), 
		this->parent_, SLOT( notify_picked_up_group_size( int ) ) );

	connect( new_group, SIGNAL( picked_up_layer_size_signal_( int ) ), 
		this->parent_, SLOT( notify_groups_of_picked_up_layer_size( int ) ) );

	return new_group;
}

void LayerManagerWidgetPrivate::reset()
{
	GroupWidgetMap::iterator it = this->group_map_.begin();
	while ( it != this->group_map_.end() )
	{
		LayerGroupWidget* widget = it->second;
		this->group_layout_->removeWidget( widget );
		widget->deleteLater();
		++it;
	}
	this->group_map_.clear();
}

void LayerManagerWidgetPrivate::HandleLayerInserted( qpointer_type qpointer, 
													LayerHandle layer, bool new_group )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind( 
		&LayerManagerWidgetPrivate::handle_layer_inserted, qpointer.data(), layer, new_group ) ) );
}

void LayerManagerWidgetPrivate::HandleLayersDeleted( qpointer_type qpointer, 
	std::vector< std::string > affected_groups, bool groups_deleted )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind( 
		&LayerManagerWidgetPrivate::handle_layers_deleted, qpointer.data(),
		affected_groups, groups_deleted ) ) );
}

void LayerManagerWidgetPrivate::HandleLayersReordered( qpointer_type qpointer, 
													  std::string group_id )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind( 
		&LayerManagerWidgetPrivate::handle_layers_reordered, qpointer.data(), group_id ) ) );
}

void LayerManagerWidgetPrivate::HandleGroupsReordered( qpointer_type qpointer )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind( 
		&LayerManagerWidgetPrivate::handle_groups_reordered, qpointer.data() ) ) );
}

void LayerManagerWidgetPrivate::PreLoadStates( qpointer_type qpointer )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind(
		&LayerManagerWidgetPrivate::pre_load_states, qpointer.data() ) ) );
}

void LayerManagerWidgetPrivate::PostLoadStates( qpointer_type qpointer )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind(
		&LayerManagerWidgetPrivate::post_load_states, qpointer.data() ) ) );
}

void LayerManagerWidgetPrivate::HandleReset( qpointer_type qpointer )
{
	Core::Interface::PostEvent( QtUtils::CheckQtPointer( qpointer, boost::bind(
		&LayerManagerWidgetPrivate::reset, qpointer.data() ) ) );
}

//////////////////////////////////////////////////////////////////////////
// Class LayerManagerWidget
//////////////////////////////////////////////////////////////////////////

LayerManagerWidget::LayerManagerWidget( QWidget* parent ) :
	QScrollArea( parent )
{
	this->private_ = new LayerManagerWidgetPrivate( this );

	// Customize the settings for the scroll area
	this->setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
	this->setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded );
	this->setContentsMargins( 1,1,1,1 );
	this->setWidgetResizable( true );
	
	this->setFrameShape( QFrame::NoFrame );
	this->setFrameShadow( QFrame::Plain );
	this->setLineWidth( 0 );

	// Define a new central widget inside the scroll area
	this->private_->main_ = new QWidget( parent );
	this->setWidget( this->private_->main_ );
	
	// Setup the spacing between the groups
	this->private_->group_layout_ = new QVBoxLayout( this->private_->main_ );
	this->private_->group_layout_->setSpacing( 2 );
	this->private_->group_layout_->setContentsMargins( 1,1,1,1 );
	this->private_->group_layout_->setAlignment( Qt::AlignTop );
	
	this->private_->main_->setSizePolicy( QSizePolicy::MinimumExpanding, QSizePolicy::Preferred );
	this->private_->main_->setAcceptDrops( true );

	LayerManagerWidgetPrivate::qpointer_type qpointer( this->private_ );

	{
		// Connect the signals from the LayerManager to the GUI

		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

		this->add_connection( LayerManager::Instance()->layer_inserted_signal_.connect( 
			boost::bind( &LayerManagerWidgetPrivate::HandleLayerInserted, qpointer, _1, _2 ) ) );
		this->add_connection( LayerManager::Instance()->layers_deleted_signal_.connect( 
			boost::bind( &LayerManagerWidgetPrivate::HandleLayersDeleted, qpointer, _2, _3 ) ) );
		this->add_connection( LayerManager::Instance()->layers_reordered_signal_.connect(
			boost::bind( &LayerManagerWidgetPrivate::HandleLayersReordered, qpointer, _1 ) ) );
		this->add_connection( LayerManager::Instance()->groups_reordered_signal_.connect(
			boost::bind( &LayerManagerWidgetPrivate::HandleGroupsReordered, qpointer ) ) );
		this->add_connection( Core::StateEngine::Instance()->pre_load_states_signal_.connect( 
			boost::bind( &LayerManagerWidgetPrivate::PreLoadStates, qpointer ) ) );
		this->add_connection( Core::StateEngine::Instance()->post_load_states_signal_.connect( 
			boost::bind( &LayerManagerWidgetPrivate::PostLoadStates, qpointer ) ) );
		this->add_connection( Core::Application::Instance()->reset_signal_.connect(
			boost::bind( &LayerManagerWidgetPrivate::HandleReset, qpointer ) ) );
	}
	
	// Add any layers that may have been added before the GUI was initialized
	this->private_->update_group_widgets();
}

LayerManagerWidget::~LayerManagerWidget()
{
	this->disconnect_all();
}
	
// Drag and Drop Functions	
void LayerManagerWidget::prep_layers_for_drag_and_drop( bool move_time )
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	for( GroupWidgetMap::iterator it = this->private_->group_map_.begin(); 
		it != this->private_->group_map_.end(); ++it )
	{
		LayerGroupHandle group = LayerManager::Instance()->find_group( ( *it ).first );
		if( group && group->group_widget_expanded_state_->get() )
		{	
			( *it ).second->prep_layers_for_drag_and_drop( move_time );
		}
	}
}

void LayerManagerWidget::prep_groups_for_drag_and_drop( bool move_time )
{
	for( GroupWidgetMap::iterator it = this->private_->group_map_.begin(); 
		it != this->private_->group_map_.end(); ++it )
	{
		( *it ).second->prep_for_animation( move_time );
	}
}
	
void LayerManagerWidget::notify_picked_up_group_size( int group_size )
{
	for( GroupWidgetMap::iterator it = this->private_->group_map_.begin(); 
		it != this->private_->group_map_.end(); ++it )
	{
		( *it ).second->set_picked_up_group_size( group_size );
	}
}

void LayerManagerWidget::notify_groups_of_picked_up_layer_size( int layer_size )
{
	for( GroupWidgetMap::iterator it = this->private_->group_map_.begin(); 
		it != this->private_->group_map_.end(); ++it )
	{
		( *it ).second->notify_picked_up_layer_size( layer_size );
	}
}

}  // end namespace Seg3D
