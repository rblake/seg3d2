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

#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/LayerManager/LayerAction.h>

namespace Seg3D
{

class LayerActionPrivate
{
public:
	// VALIDATE_LAYER_DEPENDENDIES:
	// This function serves two purposes:
	// (1) It checks if any layer ids are provenanceids. If this is the case all the provenance ids
	// are translated into normal layer ids. 
	// (2) It checks whether all the layerids or the provenanceids exist in the current layer
	// layermanager
	// NOTE: This function should be run on at the start of the validate phase of the layer action
	bool validate_layer_dependencies( ActionContextHandle& context );


};

bool LayerActionPrivate::validate_layer_dependencies( ActionContextHandle& context )
{
	

}

LayerAction::LayerAction() :
	private_( new LayerActionPrivate )

bool LayerAction::validate( ActionContextHandle& context )
{
	this->private
}




} // end namespace Seg3D
