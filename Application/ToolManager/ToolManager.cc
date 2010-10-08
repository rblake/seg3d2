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

#include <tinyxml.h>

#include <Core/Application/Application.h>
#include <Core/Interface/Interface.h>
#include <Core/State/StateIO.h>

#include <Application/Session/Session.h>
#include <Application/Tool/ToolFactory.h>
#include <Application/ToolManager/ToolManager.h>

#include <Application/ToolManager/Actions/ActionOpenTool.h>
#include <Application/ToolManager/Actions/ActionCloseTool.h>
#include <Application/ToolManager/Actions/ActionActivateTool.h>

#include <Application/ViewerManager/ViewerManager.h>

namespace Seg3D
{

CORE_SINGLETON_IMPLEMENTATION( ToolManager );

class ToolManagerPrivate
{
public:
	bool handle_mouse_enter( ViewerHandle viewer, int x, int y );
	bool handle_mouse_leave( ViewerHandle viewer );
	bool handle_mouse_move( ViewerHandle viewer, const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );
	bool handle_mouse_press( ViewerHandle viewer, const Core::MouseHistory& mouse_history, 
		int button, int buttons, int modifiers );
	bool handle_mouse_release( ViewerHandle viewer, const Core::MouseHistory& mouse_history,
		int button, int buttons, int modifiers );
	bool handle_wheel( ViewerHandle viewer, int delta, int x, int y, int buttons, int modifiers );
	bool handle_key_press( ViewerHandle viewer, int key, int modifiers );

	void update_viewers( bool redraw_2d, bool redraw_3d );

	// All the open tools are stored in this hash map
	ToolManager::tool_list_type tool_list_;
	ToolHandle active_tool_;
	ViewerHandle focus_viewer_;

	Core::StateStringHandle active_tool_state_;
};

bool ToolManagerPrivate::handle_mouse_enter( ViewerHandle viewer, int x, int y )
{
	this->focus_viewer_ = viewer;

	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_mouse_enter( viewer, x, y );
	}
	return false;
}

bool ToolManagerPrivate::handle_mouse_leave( ViewerHandle viewer )
{
	this->focus_viewer_.reset();

	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_mouse_leave( viewer );
	}
	return false;
}

bool ToolManagerPrivate::handle_mouse_move( ViewerHandle viewer, 
										   const Core::MouseHistory& mouse_history, 
										   int button, int buttons, int modifiers )
{
	// If there is no focus viewer, simulate a mouse enter event.
	// NOTE: It is not safe to pass this mouse move event to the active tool, because
	// the mouse history information might not be correct.
	if ( !this->focus_viewer_ )
	{
		return this->handle_mouse_enter( viewer, mouse_history.current_.x_, 
			mouse_history.current_.y_ );
	}
	
	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_mouse_move( viewer, mouse_history, 
			button, buttons, modifiers );
	}
	return false;
}

bool ToolManagerPrivate::handle_mouse_press( ViewerHandle viewer, 
											const Core::MouseHistory& mouse_history, 
											int button, int buttons, int modifiers )
{
	this->focus_viewer_ = viewer;

	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_mouse_press( viewer, mouse_history, 
			button, buttons, modifiers );
	}
	return false;
}

bool ToolManagerPrivate::handle_mouse_release( ViewerHandle viewer,
											  const Core::MouseHistory& mouse_history, 
											  int button, int buttons, int modifiers )
{
	this->focus_viewer_ = viewer;

	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_mouse_release( viewer, mouse_history, 
			button, buttons, modifiers );
	}
	return false;
}

bool ToolManagerPrivate::handle_wheel( ViewerHandle viewer, int delta, 
									  int x, int y, int buttons, int modifiers )
{
	this->focus_viewer_ = viewer;

	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_wheel( viewer, delta, x, y, buttons, modifiers );
	}
	return false;
}

bool ToolManagerPrivate::handle_key_press( ViewerHandle viewer, int key, int modifiers )
{
	this->focus_viewer_ = viewer;

	ToolHandle active_tool;
	{
		ToolManager::lock_type lock( ToolManager::Instance()->get_mutex() );
		active_tool = this->active_tool_;
	}
	if ( active_tool )
	{
		return active_tool->handle_key_press( viewer, key, modifiers );
	}
	return false;
}

void ToolManagerPrivate::update_viewers( bool redraw_2d, bool redraw_3d )
{
	size_t num_of_viewers = ViewerManager::Instance()->number_of_viewers();
	for ( size_t i = 0; i < num_of_viewers; i++ )
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_viewer( i );
		if ( viewer->is_volume_view() )
		{
			if ( redraw_3d )
			{
				viewer->redraw();
			}
		}
		else if ( redraw_2d )
		{
			viewer->redraw_overlay();
		}
	}
}

ToolManager::ToolManager() :
	StateHandler( "toolmanager", false ),
	private_( new ToolManagerPrivate )
{
	this->add_state( "active_tool", this->private_->active_tool_state_, Tool::NONE_OPTION_C );

	// Register mouse event handlers for all the viewers
	size_t num_of_viewers = ViewerManager::Instance()->number_of_viewers();
	for ( size_t i = 0; i < num_of_viewers; i++ )
	{
		ViewerHandle viewer = ViewerManager::Instance()->get_viewer( i );
		viewer->set_mouse_enter_handler( boost::bind( &ToolManagerPrivate::handle_mouse_enter,
			this->private_, _1, _2, _3 ) );
		viewer->set_mouse_leave_handler( boost::bind( &ToolManagerPrivate::handle_mouse_leave,
			this->private_, _1 ) );
		viewer->set_mouse_move_handler( boost::bind( &ToolManagerPrivate::handle_mouse_move,
			this->private_, _1, _2, _3, _4, _5 ) );
		viewer->set_mouse_press_handler( boost::bind( &ToolManagerPrivate::handle_mouse_press,
			this->private_, _1, _2, _3, _4, _5 ) );
		viewer->set_mouse_release_handler( boost::bind( &ToolManagerPrivate::handle_mouse_release,
			this->private_, _1, _2, _3, _4, _5 ) );
		viewer->set_wheel_event_handler( boost::bind( &ToolManagerPrivate::handle_wheel,
			this->private_, _1, _2, _3, _4, _5, _6 ) );
		viewer->set_key_press_event_handler( boost::bind( &ToolManagerPrivate::handle_key_press,
			this->private_, _1, _2, _3 ) );
	}
}

ToolManager::~ToolManager()
{
	this->disconnect_all();
}

// THREAD-SAFETY:
// Only ActionOpenTool calls this function and this action is only run on the
// application thread. Hence the function is always executed by the same thread.

bool ToolManager::open_tool( const std::string& tool_type, std::string& new_toolid )
{
	// Step (1): Make the function thread safe
	lock_type lock( this->get_mutex() );

	// Step (2): Add an entry in the debug log
	CORE_LOG_MESSAGE( std::string( "Open tool: " ) + tool_type );

	// Step (4): Build the tool using the factory. This will generate the default
	// settings.
	ToolHandle tool;

	if ( !( ToolFactory::Instance()->create_tool( tool_type, tool ) ) )
	{
		CORE_LOG_ERROR( std::string( "Could not create tool of type: '" ) + tool_type + "'" );
		return false;
	}

	// Step (5): Add the tool id to the tool and add the tool to the list
	new_toolid = tool->toolid();
	this->private_->tool_list_[ new_toolid ] = tool;
	
	// Step (6): Signal any observers (UIs) that the tool has been opened
	open_tool_signal_( tool );

	// Set the tool to be active
	this->activate_tool( new_toolid );

	// All done
	return true;
}

// THREAD-SAFETY:
// Only ActionCloseTool calls this function and this action is only run on the
// application thread. Hence the function is always executed by the same thread.

void ToolManager::close_tool( const std::string& toolid )
{
	// Step (1): Make the function thread safe
	lock_type lock( this->get_mutex() );

	// Step (2): Add an entry in the debug log
	CORE_LOG_MESSAGE( std::string( "Close tool: " ) + toolid );

	// Step (3): Find the tool in the list.
	tool_list_type::iterator it = this->private_->tool_list_.find( toolid );
	if ( it == this->private_->tool_list_.end() )
	{
		CORE_LOG_ERROR( std::string( "Toolid '" + toolid + "' does not exist" ) );
		return;
	}

	// Step (4): Get the tool from the iterator
	ToolHandle tool = ( *it ).second;
	if ( tool == this->private_->active_tool_ )
	{
		// Call the tool deactivate function that will unregister the current
		// tool bindings
		tool->deactivate();

		// Set no tool as active
		this->private_->active_tool_.reset();
	}

	// Step (5): Move the tool from the list. The tool handle still persists
	// and will be removed after the signal has been posted.
	this->private_->tool_list_.erase( it );

	// Step (6): Run the function in the tool that cleans up the parts that
	// need to be cleaned up on the interface thread.
	tool->close();

	// Step (7): Signal that the tool will be closed.
	close_tool_signal_( tool );

	if ( !this->private_->active_tool_ )
	{
		if ( !this->private_->tool_list_.empty() )
		{
			this->private_->active_tool_ = ( *this->private_->tool_list_.begin() ).second;
			this->private_->active_tool_->activate();
			this->activate_tool_signal_( this->private_->active_tool_ );
		}
		else if ( tool->has_2d_visual() || tool->has_3d_visual() )
		{
			this->private_->update_viewers( tool->has_2d_visual(), tool->has_3d_visual() );
		}
	}
}

// THREAD-SAFETY:
// Only ActionActivateTool calls this function and this action is only run on the
// application thread. Hence the function is always executed by the same thread.

void ToolManager::activate_tool( const std::string& toolid )
{
	// Step (1): Make the function thread safe
	lock_type lock( this->get_mutex() );

	// Step (2): Add an entry in the debug log
	CORE_LOG_DEBUG( std::string( "Activate tool: " ) + toolid );

	// Step (3): Check if anything needs to be done
	if ( this->private_->active_tool_ && 
		this->private_->active_tool_->toolid() == toolid )
	{
		return;
	}

	// Step (4): Find new active tool
	tool_list_type::iterator it = this->private_->tool_list_.find( toolid );
	if ( it == this->private_->tool_list_.end() )
	{
		return;
	}

	// Step (4): Deactivate the current active tool if it exists, and activate the new one
	ToolHandle old_tool = this->private_->active_tool_;
	this->private_->active_tool_ = ( *it ).second;
	if ( old_tool )
	{
		old_tool->deactivate();
	}	
	this->private_->active_tool_->activate();

	// Step (5): Update viewers if necessary.
	bool redraw_2d = this->private_->active_tool_->has_2d_visual() ||
		( old_tool && old_tool->has_2d_visual() );
	bool redraw_3d = this->private_->active_tool_->has_3d_visual() ||
		( old_tool && old_tool->has_3d_visual() );

	if ( redraw_2d || redraw_3d )
	{
		this->private_->update_viewers( redraw_2d, redraw_3d );
	}
	
	// Step (6): signal for interface
	activate_tool_signal_( ( *it ).second );
}

ToolManager::tool_list_type ToolManager::tool_list()
{
	lock_type lock( this->get_mutex() );
	return this->private_->tool_list_;
}

std::string ToolManager::active_toolid()
{
	lock_type lock( this->get_mutex() );
	if ( this->private_->active_tool_ )
	{
		return this->private_->active_tool_->toolid();
	}
	return Tool::NONE_OPTION_C;
}

ToolHandle ToolManager::get_active_tool()
{
	lock_type lock( this->get_mutex() );
	return this->private_->active_tool_;
}

bool ToolManager::pre_save_states( Core::StateIO& state_io )
{
	this->private_->active_tool_state_->set( this->active_toolid() );
	return true;
}

bool ToolManager::post_save_states( Core::StateIO& state_io )
{
	TiXmlElement* tm_element = state_io.get_current_element();
	assert( this->get_statehandler_id() == tm_element->Value() );
	TiXmlElement* tools_element = new TiXmlElement( "tools" );
	tm_element->LinkEndChild( tools_element );

	state_io.push_current_element();
	state_io.set_current_element( tools_element );

	tool_list_type::iterator it = this->private_->tool_list_.begin();
	tool_list_type::iterator it_end = this->private_->tool_list_.end();
	while ( it != it_end )
	{
		( *it ).second->save_states( state_io );
		++it;
	}

	state_io.pop_current_element();

	return true;
}

bool ToolManager::post_load_states( const Core::StateIO& state_io )
{
	const TiXmlElement* tools_element = state_io.get_current_element()->
		FirstChildElement( "tools" );
	if ( tools_element == 0 )
	{
		return false;
	}
	
	state_io.push_current_element();
	state_io.set_current_element( tools_element );

	bool success = true;
	const TiXmlElement* tool_element = tools_element->FirstChildElement();
	while ( tool_element != 0 )
	{
		std::string toolid( tool_element->Value() );
		if ( this->open_tool( toolid, toolid ) )
		{
			ToolHandle tool = this->get_tool( toolid );
			success &= tool->load_states( state_io );
		}
		tool_element = tool_element->NextSiblingElement();
	}

	state_io.pop_current_element();

	if( this->private_->active_tool_state_->get() != Tool::NONE_OPTION_C )
	{
		this->activate_tool( this->private_->active_tool_state_->get() );
	}
	else if ( !this->private_->tool_list_.empty() )
	{
		this->activate_tool( ( *this->private_->tool_list_.begin() ).first );
	}

	return true;
}

bool ToolManager::pre_load_states( const Core::StateIO& state_io )
{
	return this->delete_all();
}

ToolHandle ToolManager::get_tool( const std::string& toolid )
{
	tool_list_type::iterator it = this->private_->tool_list_.find( toolid );
	if ( it != this->private_->tool_list_.end() )
	{
		return ( *it ).second;
	}
	
	return ToolHandle();
}

void ToolManager::get_tool_names( std::vector< ToolIDNamePair >& tool_names )
{
	lock_type lock( this->get_mutex() );
	tool_list_type::iterator it = this->private_->tool_list_.begin();
	tool_list_type::iterator it_end = this->private_->tool_list_.end();
	while ( it != it_end )
	{
		ToolHandle tool = ( *it ).second;
		tool_names.push_back( std::make_pair( tool->toolid(), tool->tool_name() ) );
		it++;
	}
}

bool ToolManager::delete_all()
{
	tool_list_type::iterator it = this->private_->tool_list_.begin();
	while( it != this->private_->tool_list_.end() )
	{
		( *it ).second->invalidate();
		tool_list_type::iterator temp_it = it;
		it++;
		this->close_tool( ( *temp_it ).first );
	}
	return true;
}

int ToolManager::get_session_priority()
{
	return SessionPriority::TOOL_MANAGER_PRIORITY_E;
}

} // end namespace Seg3D
