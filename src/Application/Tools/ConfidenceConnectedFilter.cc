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
#include <Application/Filters/Actions/ActionConfidenceConnectedFilter.h>
#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/ConfidenceConnectedFilter.h>
#include  <Application/ViewerManager/ViewerManager.h>

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ConfidenceConnectedFilter )

namespace Seg3D
{

ConfidenceConnectedFilter::ConfidenceConnectedFilter( const std::string& toolid ) :
	SeedPointsTool( Core::VolumeType::DATA_E, toolid )
{
	this->add_state( "iterations", this->iterations_state_, 3, 1, 100, 1 );
	this->add_state( "multiplier", this->multiplier_state_, 2.0, 0.1, 100.0, 0.1 );
}

ConfidenceConnectedFilter::~ConfidenceConnectedFilter()
{
	this->disconnect_all();
}

void ConfidenceConnectedFilter::execute( Core::ActionContextHandle context )
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
	ActionConfidenceConnectedFilter::Dispatch( 
		context, this->target_layer_state_->get(),
		this->seed_points_state_->get(), 
		static_cast< unsigned int >( this->iterations_state_->get() ), 
		this->multiplier_state_->get() );
}

} // end namespace Seg3D
