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

// Application includes
#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/PaintTool.h>
#include <Application/LayerManager/LayerManager.h>

namespace Seg3D
{

// Register the tool into the tool factory
SCI_REGISTER_TOOL(PaintTool)

PaintTool::PaintTool( const std::string& toolid ) :
	Tool( toolid )
{
	// Need to set ranges and default values for all parameters
	add_state( "target", this->target_layer_state_, "<none>" );
	add_state( "mask", this->mask_layer_state_, "<none>" );
	add_state( "brush_radius", this->brush_radius_state_, 23, 1, 250, 1 );
	add_state( "upper_threshold", this->upper_threshold_state_, 1.0, 00.0, 1.0, 0.01 );
	add_state( "lower_threshold", this->lower_threshold_state_, 0.0, 00.0, 1.0, 0.01 );
	add_state( "erase", this->erase_state_, false );
	
	this->handle_layers_changed();

	// Add constaints, so that when the state changes the right ranges of
	// parameters are selected
	this->add_connection ( this->target_layer_state_->value_changed_signal_.connect( boost::bind( &PaintTool::target_constraint,
	    this, _1 ) ) );
	this->add_connection ( this->mask_layer_state_->value_changed_signal_.connect( boost::bind( &PaintTool::mask_constraint,
	    this, _1 ) ) );
	
	this->add_connection ( LayerManager::Instance()->layers_changed_signal_.connect(
		boost::bind( &PaintTool::handle_layers_changed, this ) ) );

}

PaintTool::~PaintTool()
{
	disconnect_all();
}

void PaintTool::handle_layers_changed()
{
	std::vector< LayerHandle > target_layers;
	LayerManager::Instance()->get_layers( target_layers );
	bool target_found = false;
	bool mask_found = false;
	
	for( int i = 0; i < static_cast< int >( target_layers.size() ); ++i )
	{
		if( ( this->target_layer_state_->get() == "<none>" ) && ( target_layers[i]->type() == 
																 Core::VolumeType::DATA_E ) )
		{
			this->target_layer_state_->set( target_layers[i]->get_layer_name(), 
				Core::ActionSource::NONE_E );
		}
		
		if( ( this->mask_layer_state_->get() == "<none>" ) && ( target_layers[i]->type() == 
																 Core::VolumeType::MASK_E ) )
		{
			this->mask_layer_state_->set( target_layers[i]->get_layer_name(), 
				Core::ActionSource::NONE_E );
		}
		
		if( target_layers[i]->get_layer_name() == this->target_layer_state_->get() ) 
			target_found = true;
	
		if( target_layers[i]->get_layer_name() == this->mask_layer_state_->get() )
			mask_found = true;
	}
		
	if( !target_found )
		this->target_layer_state_->set( "", Core::ActionSource::NONE_E );
	
	if( !mask_found )
		this->mask_layer_state_->set( "", Core::ActionSource::NONE_E );

}

void PaintTool::target_constraint( std::string layerid )
{
}

void PaintTool::mask_constraint( std::string layerid )
{
}

void PaintTool::activate()
{
}

void PaintTool::deactivate()
{
}

} // end namespace Seg3D
