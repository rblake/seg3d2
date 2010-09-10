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
#include <Application/Filters/LayerResampler.h>
#include <Application/Filters/Actions/ActionResample.h>
#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerGroup.h>

namespace Seg3D
{

class LayerResamplerPrivate
{
public:
	void handle_kernel_changed( std::string kernel_name );

	LayerResampler* tool_;
	LayerHandle src_layer_;
	LayerGroupHandle dst_group_;
};

void LayerResamplerPrivate::handle_kernel_changed( std::string kernel_name )
{
	this->tool_->has_params_state_->set( kernel_name == ActionResample::GAUSSIAN_C );
}

//////////////////////////////////////////////////////////////////////////
// Class LayerResampler
//////////////////////////////////////////////////////////////////////////

LayerResampler::LayerResampler( LayerHandle src_layer, LayerGroupHandle dst_group ) :
	Core::StateHandler( "layerresampler", false ),
	private_( new LayerResamplerPrivate )
{
	this->private_->tool_ = this;
	this->private_->src_layer_ = src_layer;
	this->private_->dst_group_ = dst_group;

	std::vector< Core::OptionLabelPair > padding_values;
	padding_values.push_back( std::make_pair( ActionResample::ZERO_C, "0" ) );
	padding_values.push_back( std::make_pair( ActionResample::MIN_C, "Minimum Value" ) );
	padding_values.push_back( std::make_pair( ActionResample::MAX_C, "Maximum Value" ) );
	this->add_state( "pad_value", this->padding_value_state_, ActionResample::ZERO_C, 
		padding_values );

	std::vector< Core::OptionLabelPair > kernels;
	kernels.push_back( std::make_pair( ActionResample::BOX_C, "Box" ) );
	kernels.push_back( std::make_pair( ActionResample::TENT_C, "Tent" ) );
	kernels.push_back( std::make_pair( ActionResample::CUBIC_CR_C, "Cubic (Catmull-Rom)" ) );
	kernels.push_back( std::make_pair( ActionResample::CUBIC_BS_C, "Cubic (B-spline)" ) );
	kernels.push_back( std::make_pair( ActionResample::QUARTIC_C, "Quartic" ) );
	kernels.push_back( std::make_pair( ActionResample::GAUSSIAN_C, "Gaussian" ) );
	this->add_state( "kernel", this->kernel_state_, ActionResample::BOX_C, kernels );

	this->add_state( "sigma", this->gauss_sigma_state_, 1.0, 1.0, 100.0, 0.01 );
	this->add_state( "cutoff", this->gauss_cutoff_state_, 1.0, 1.0, 100.0, 0.01 );
	this->add_state( "has_params", this->has_params_state_, false );

	this->add_connection( this->kernel_state_->value_changed_signal_.connect(
		boost::bind( &LayerResamplerPrivate::handle_kernel_changed, this->private_, _2 ) ) );
}

LayerResampler::~LayerResampler()
{
	this->disconnect_all();
}

void LayerResampler::execute( Core::ActionContextHandle context )
{
	Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );

	std::vector< std::string > layer_ids( 1, this->private_->src_layer_->get_layer_id() );
	ActionResample::Dispatch( context, layer_ids,
		this->private_->dst_group_->get_grid_transform(), 
		this->padding_value_state_->get(),
		this->kernel_state_->get(), this->gauss_sigma_state_->get(),
		this->gauss_cutoff_state_->get(), true );
}

} // end namespace Seg3D
