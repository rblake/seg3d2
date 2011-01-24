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
 
#ifndef APPLICATION_LAYERMANAGER_LAYERACTION_H
#define APPLICATION_LAYERMANAGER_LAYERACTION_H 
 
#include <Core/Action/Action.h> 
#include <Application/Layer/LayerFWD.h> 

namespace Seg3D
{

// Forward declarations
class LayerAction;
class LayerActionPrivate;
typedef boost::shared_ptr<LayerAction> LayerActionHandle;
typedef boost::shared_ptr<LayerActionPrivate> LayerActionPrivateHandle;

class LayerAction : public Core::Action
{
	// -- constructor / destructor
public:
	// Constructor that generates the private class
	LayerAction();

	// Virtual destructor for memory management of derived classes
	virtual ~LayerAction();

public:
	// VALIDATE_LAYER_ACTION:
	// Each action needs to be validated just before it is posted. This way we
	// enforce that every action that hits the main post_action signal will be
	// a valid action to execute.
	virtual bool validate_layer_action( ActionContextHandle& context ) = 0;

	// RUN_LAYER_ACTION:
	// Each action needs to have this piece implemented. It spells out how the
	// action is run. It returns whether the action was successful or not.
	// NOTE: In case of an asynchronous action, the return value is ignored and the
	// program relies on report_done() from the context to be triggered when
	// the asynchronous part has finished. In any other case the ActionDispatcher
	// will issue the report_done() when run returns.
	virtual bool run_layer_action( ActionContextHandle& context, ActionResultHandle& result ) = 0;

public:
	// NOTE: These functions are implemented in the LayerAction and should *NOT* be overloaded in
	// the derived class. The functions validate_layer_action and run_layer_action should be called
	// instead.
	
	virtual bool validate( ActionContextHandle& context );
	virtual bool run( ActionContextHandle& context, ActionResultHandle& result );

	// -- deal with dependencies for provenance --
protected:
	// NOTE: Dependency lists need to inserted into a ProvenanceRecord and hence this function can
	// only be called from within run_layer_action or validate_layer_action
	
	// GET_INPUT_PROVENANCE_IDS:
	// Get the dependencies that this action is depending on
	ProvenanceIDList get_input_provenance_ids();
	
	// -- deal with new provenance ids for output --
public: 
	// SET_OUTPUT_PROVENANCE_IDS:
	// Set the output provenance ids
	void set_output_provenance_ids( const ProvenanceIDList& provenance_ids );
	
	// GENERATE_OUTPUT_PROVENANCE_ID:
	// Get the provenance id of output layer indexed by index. If no provenance id was assigned
	// a new one is created
	ProvenanceID generate_output_provenance_id( size_t index = 0 );

	// GET_OUTPUT_PROVENANCE_IDS:
	// Get all the assigned provenance ids
	ProvenanceIDList get_output_provenance_ids();
	
	// -- internals --
private:	
	LayerActionPrivateHandle private_;
};

} // end namespace Seg3D

#define LAYER_ACTION(definition_string) \
CORE_ACTION_INTERNAL( definition_string CORE_ACTION_CHANGES_PROJECT_DATA() \
CORE_ACTION_CHANGES_PROVENANCE_DATA() CORE_ACTION_IS_UNDOABLE(), Core::LayerActionInfo )

#define LAYER_ACTION_TYPE( name, description ) \
CORE_ACTION_TYPE( name, description )

#define LAYER_ACTION_ARGUMENT( name, description ) \
CORE_ACTION_ARGUMENT( name, description ) 

#define LAYER_ACTION_LAYERID_ARGUMENT( name, description ) \
CORE_ACTION_ARGUMENT( name, description ) \
CORE_ACTION_PROVENANCE_ID_PARAMETER( name )


#endif
