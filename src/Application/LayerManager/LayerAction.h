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
	// -- constructor --
public:
	// Constructor that generates the private class
	LayerAction();

	// -- functionality for setting parameter list --
public:
	// ADD_PARAMETER:
	// Add a parameter to the internal database
	template< class T >
	void add_parameter( const T& parameter )
	{
		this->add_parameter_internal(
			Core::ActionParameterBaseHandle( new Core::ActionParameter<T>( &parameter ) ) );
	}
	
	// ADD_PARAMETER:
	// Specialized function for adding and registering layer id input parameters
	void add_parameter( const InputLayerID& layer_id );
	
	// ADD_PARAMETER:
	// Specialized function for adding and registering layer ids input parameters
	void add_parameter( const std::vector<InputLayerID>& layer_ids );

	// -- translate provenance information --
public:
	// TRANSLATE:
	// Some actions need to be translated before they can be validated. Translate takes
	// care of most provenance related issue, by for example translating the provenance
	// information into real action information. This function is called before validate
	// NOTE: This function is *not* const and may alter the values of the parameters
	//       and correct faulty input.
	virtual bool translate( ActionContextHandle& context );
	
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

#endif
