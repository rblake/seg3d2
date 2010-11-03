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

#ifndef APPLICATION_LAYERMANAGER_LAYERUNDOBUFFER_H
#define APPLICATION_LAYERMANAGER_LAYERUNDOBUFFER_H

// Application includes
#include <Application/LayerManager/LayerUndoBufferItem.h>

namespace Seg3D
{

// LayerUndoBuffer
// This singleton class keeps a list of the actions that can be undone or redone.

// Internals for the LayerUndoBuffer singleton
class LayerUndoBufferPrivate;
typedef boost::shared_ptr<LayerUndoBufferPrivate> LayerUndoBufferPrivateHandle;

class LayerUndoBuffer : public Core::ConnectionHandler
{
	CORE_SINGLETON( LayerUndoBuffer );

	// -- Constructor/Destructor --
private:
	LayerUndoBuffer();
	virtual ~LayerUndoBuffer();

	// -- undo interface --
public:

	// INSERT_UNDO_ITEM:
	// Insert a new undo item in the queue
	// NOTE: The action context is needed to verify whether it is inserted from the undo buffer
	// itself or whether the undo item was created in a normal action.
	void insert_undo_item( Core::ActionContextHandle context, 
		LayerUndoBufferItemHandle undo_item );

	// UNDO:
	// Undo the top item of the stack
	bool undo( Core::ActionContextHandle context  ); 
	
	// REDO:
	// redo the top item of the stack
	bool redo( Core::ActionContextHandle context );
	
	// RESET_UNDO_BUFFER:
	// Reset the buffer to its initial setting
	void reset_undo_buffer();
	
	// GET_UNDO_TAG:
	// Get the tag from the action stored on top of the undo stack.
	// NOTE: if an empty string is returned the undo stack is empty.
	std::string get_undo_tag() const;
	
	// GET_REDO_TAG:
	// Get the tag from the action stored on top of the redo stack.
	// NOTE: if an empty string is returned the redo stack is empty.
	std::string get_redo_tag() const;

	// HAS_UNDO:
	// Check whether there is something to undo
	bool has_undo() const;
	
	// HAS_REDO:
	// Check whether there is something to redo
	bool has_redo() const;

	// -- signals --
public:
	typedef boost::signals2::signal< void ( std::string ) > update_undo_tag_signal_type;
	typedef boost::signals2::signal< void ( std::string ) > update_redo_tag_signal_type;

	// UPDATE_UNDO_TAG_SIGNAL:
	// This signal is triggered when a new undo item is on top of the undo stack
	update_undo_tag_signal_type update_undo_tag_signal_;
	
	// UPDATE_REDO_TAG_SIGNAL:
	// This signal is triggered when a new redo item is on top of the redo stack
	update_redo_tag_signal_type update_redo_tag_signal_;
	
	// -- internals --
private:
	// Handle to internals
	LayerUndoBufferPrivateHandle private_;
};

} // end namespace Seg3D

#endif
