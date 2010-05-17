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

// TinyXML includes
#include <Externals/tinyxml/tinyxml.h>

#include <Core/State/StateHandler.h>
#include <Core/State/StateEngine.h>
#include <Core/Application/Application.h>

namespace Core
{

typedef std::map< std::string, StateBaseHandle > state_map_type;

class StateHandlerPrivate
{
public:
	// The id of this state handler
	std::string statehandler_id_;

	// The database with the actual states 
	state_map_type state_map_;
};


StateHandler::StateHandler( const std::string& type_str, bool auto_id )
{
	this->private_ = new StateHandlerPrivate;
	this->private_->statehandler_id_ = StateEngine::Instance()->
		register_state_handler( type_str, this, auto_id );

}

StateHandler::~StateHandler()
{
	this->disconnect_all();
	StateEngine::Instance()->remove_state_handler( this->private_->statehandler_id_ );
	delete this->private_;
}
		

bool StateHandler::add_statebase( StateBaseHandle state )
{
	// Step (1): Get unique state id
	std::string stateid = state->stateid();

	// Step (2): Link with statehandler
	this->add_connection( state->state_changed_signal_.connect( boost::bind(
	    &StateHandler::handle_state_changed, this ) ) );

	// Step (3): Add the state to the map
	this->private_->state_map_[ stateid ] = state;

	return true;
}

void StateHandler::state_changed()
{
	// default function is to do nothing
}


const std::string& StateHandler::get_statehandler_id() const
{
	return ( this->private_->statehandler_id_ );
}

std::string StateHandler::create_state_id( const std::string& key ) const
{
	return this->get_statehandler_id() + "::" + key;
}

bool StateHandler::get_state( const std::string& state_id, StateBaseHandle& state )
{
	state_map_type::iterator it = this->private_->state_map_.find( state_id );
	if ( it != this->private_->state_map_.end() )
	{
		state = ( *it ).second;
		return true;
	}

	state.reset();
	return false;
}

bool StateHandler::get_state( const size_t idx, StateBaseHandle& state )
{
	state_map_type::iterator it = this->private_->state_map_.begin();
	
	if ( idx >= this->private_->state_map_.size() ) return false;
	
	std::advance( it, idx );
	
	state = ( *it ).second;
	return true;
}

bool StateHandler::save_states( boost::filesystem::path path )
{
	state_map_type::iterator it = this->private_->state_map_.begin();
	state_map_type::iterator it_end = this->private_->state_map_.end();

	// XML declaration and version number
	TiXmlDocument doc;  
 	TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );  
	doc.LinkEndChild( decl );  

	// Set Seg3D2 as our root
	TiXmlElement * root = new TiXmlElement( "Seg3D2" );  
	doc.LinkEndChild( root );  
 
	// We will use our statehandler_id as the parent for its state values
	TiXmlElement* preferences = new TiXmlElement( this->private_->statehandler_id_.c_str() );  
	root->LinkEndChild( preferences );  

	TiXmlElement* state_value;

	while ( it != it_end )
	{
		// Here we will add all the state values 
		state_value = new TiXmlElement( ( *it ).second->stateid().c_str() );
		state_value->LinkEndChild( new TiXmlText( ( *it ).second->export_to_string().c_str() ));  
		preferences->LinkEndChild( state_value );
		++it;
	}

	// Finally we will save our XML to the specified file
	doc.SaveFile( path.string().c_str() );
	
	return true;
}


bool StateHandler::load_states( boost::filesystem::path path )
{
	// We will load in the file from the specified path and exit if the path is invalid
	TiXmlDocument doc( path.string().c_str() );
	if ( !doc.LoadFile() ) return false;

	TiXmlHandle hDoc( &doc );
	TiXmlElement* state_values;
	TiXmlHandle hRoot(0);

	// We should have a valid root if not we will exit
	{
		state_values = hDoc.FirstChildElement().Element();
		if ( !state_values ) return false;
		hRoot = TiXmlHandle( state_values );
	}

	// Now we are expecting to get the proper statehandler_id_
	{
		state_values = hRoot.FirstChild( this->private_->
			statehandler_id_.c_str() ).FirstChild().Element();

		for( state_values; state_values; state_values = state_values->NextSiblingElement() )
		{
			// Finally we import the actual state values from the XML and import them
			std::string state_value_name = std::string( state_values->Value() );
			std::string state_value_value = std::string( state_values->GetText() );
			if( ( state_value_name != "" ) && ( state_value_value != "" ) )
				private_->state_map_[ state_value_name ]->import_from_string( state_value_value );
		}
	}
	return true;
}

void StateHandler::handle_state_changed()
{
  // Trigger the signal in the state engine
  CORE_LOG_DEBUG("Triggering state changed signal");
  StateEngine::Instance()->state_changed_signal_();
  
  // Call the local function of this state engine that handles the specifics of the derived
  // class when the state engine has changed
  state_changed();
}

size_t StateHandler::number_of_states() const
{
	return this->private_->state_map_.size();
}

} // end namespace Core
