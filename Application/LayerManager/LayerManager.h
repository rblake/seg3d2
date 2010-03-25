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

#ifndef APPLICATION_LAYERMANAGER_LAYERMANAGER_H
#define APPLICATION_LAYERMANAGER_LAYERMANAGER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 



// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

// Core includes
#include <Utils/Core/StringUtil.h>
#include <Utils/Singleton/Singleton.h>
#include <Utils/Core/Log.h>
#include <Utils/Core/Exception.h>

// Application includes
#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>
#include <Application/LayerManager/LayerScene.h>
#include <Application/State/StateHandler.h>

namespace Seg3D
{

// Forward declaration
class LayerManager;

// Class definition
class LayerManager : public StateHandler, public Utils::Singleton< LayerManager >
{

	// -- Constructor/Destructor --

private:
	friend class Utils::Singleton< LayerManager >;
	LayerManager();
	virtual ~LayerManager();

	
public:
	typedef std::list < LayerGroupHandle > group_handle_list_type;
	
	// Accessor Functions
public:
    // Functions for getting a copy of the Layers and Groups with the proper locking
	
	// TODO: get_groups() get_layers()
	void return_group_vector( std::vector< LayerGroupHandle > &vector_of_groups );
	void return_layers_vector( std::vector< LayerHandle > &vector_of_layers );
	
	LayerGroupHandle check_for_group( std::string group_id );
	
	LayerGroupHandle get_active_group();

	LayerHandle get_active_layer()
    {
        return active_layer_;
    }
	
	// Action Functions
public:
	bool insert_layer( LayerHandle layer );
	void delete_layers( LayerGroupHandle group );
	void set_active_layer( LayerHandle layer );
	
	//friend class ActionInsertLayer;

public:
	// Take an atomic snapshot of visual properties of layers for rendering in the specified viewer
	LayerSceneHandle compose_layer_scene( size_t viewer_id );
	
	// Typedef's for the Mutex
public: 
	typedef boost::recursive_mutex mutex_type;
	typedef boost::unique_lock< mutex_type > lock_type;
	
	// GET_MUTEX:
	// Get the mutex, so it can be locked by the interface that is built
	// on top of this
	mutex_type& get_mutex();

public:
	// -- Signal/Slots --
	typedef boost::signals2::signal< void( LayerHandle ) > layer_signal_type;
	typedef boost::signals2::signal< void( LayerGroupHandle ) > group_signal_type;
	typedef boost::signals2::signal< void( std::vector< LayerHandle > ) > layers_signal_type;
	// ACTIVE_LAYER_CHANGED_SIGNAL:
	// This signal is triggered after the active layer is changed
	layer_signal_type active_layer_changed_signal_; 
	
	// LAYER_INSERTED_SIGNAL:
	// This signal is triggered after a layer has been inserted
	layer_signal_type layer_inserted_signal_;
	
	// LAYER_DELETED_SIGNAL:
	// This signal is triggered after a layer has been deleted
	layer_signal_type layer_deleted_signal_;
	
	// LAYERS_DELETED_SIGNAL:
	// This signal is triggered after layers have been deleted and passes a vector of layerhandles
	layers_signal_type layers_deleted_signal_;
	
	// LAYERS_FINISHED_DELETING_SIGNAL:
	// This signal is triggered once the layers from a particular group have been deleted
	// TODO: maybe group_finished_deleting_layers_signal_
	group_signal_type layers_finished_deleting_signal_;
	
	// TODO: group_inserted_signal_ and group_deleted_signal_
	group_signal_type add_group_signal_; // Probably not needed 
	group_signal_type delete_group_signal_;
	
private:
    // list of the current groups	
	group_handle_list_type group_handle_list_;
	
	// currently active layer
	LayerHandle active_layer_;

		
};

} // end namespace seg3D

#endif

