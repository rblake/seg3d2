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

// Core includes
#include <Core/DataBlock/MaskDataBlockManager.h>
#include <Core/Parser/ArrayMathEngine.h>

// Application includes
#include <Application/Layer/LayerGroup.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/Filters/Actions/ActionArithmeticFilter.h>
#include <Application/Filters/BaseFilter.h>
#include <Application/StatusBar/StatusBar.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
// NOTE: Registration needs to be done outside of any namespace
CORE_REGISTER_ACTION( Seg3D, ArithmeticFilter )

namespace Seg3D
{

//////////////////////////////////////////////////////////////////////////
// ALGORITHM CLASS
// This class does the actual work and is run on a separate thread.
// NOTE: The separation of the algorithm into a private class is for the purpose of running the
// filter on a separate thread.
//////////////////////////////////////////////////////////////////////////

class ArithmeticFilterAlgo : public BaseFilter
{

public:
	LayerHandle dst_layer_;
	Core::ArrayMathEngine engine_;

public:

	// RUN:
	// Implementation of run of the Runnable base class, this function is called 
	// when the thread is launched.
	virtual void run();

	// GET_FITLER_NAME:
	// The name of the filter, this information is used for generating new layer labels.
	virtual std::string get_filter_name() const
	{
		return "Arithmetic";
	}
};


void ArithmeticFilterAlgo::run()
{
	std::string error;
	if ( !this->engine_.run( error ) )
	{
		StatusBar::SetMessage( Core::LogMessageType::ERROR_E, error );
		return;
	}

	if ( !this->check_abort() )
	{
		Core::DataBlockHandle output_data_block;
		this->engine_.get_data_block( ActionArithmeticFilter::RESULT_C, output_data_block );
		if ( this->dst_layer_->type() == Core::VolumeType::MASK_E )
		{
			Core::MaskDataBlockHandle mask_data_block;
			Core::MaskDataBlockManager::Convert( output_data_block, 
				this->dst_layer_->get_grid_transform(), mask_data_block );
			this->dispatch_insert_mask_volume_into_layer( this->dst_layer_, Core::MaskVolumeHandle(
				new Core::MaskVolume( this->dst_layer_->get_grid_transform(), 
				mask_data_block ) ), true );
		}
		else
		{
			this->dispatch_insert_data_volume_into_layer( this->dst_layer_, Core::DataVolumeHandle(
				new Core::DataVolume( this->dst_layer_->get_grid_transform(), output_data_block ) ), 
				false, true );
		}
	}

	this->dst_layer_->update_progress( 1.0 );
}

//////////////////////////////////////////////////////////////////////////
// Class ActionArithmeticFilterPrivate
//////////////////////////////////////////////////////////////////////////

class ActionArithmeticFilterPrivate
{
public:
	// VALIDATE:
	// Validate the filter inputs, outputs and expression prior to running filter.
	bool validate_parser( Core::ActionContextHandle& context );

	Core::ActionParameter< std::vector< std::string > > layer_ids_;
	Core::ActionParameter< std::string > expressions_;
	Core::ActionParameter< std::string > output_type_;
	Core::ActionParameter< bool > replace_;
	Core::ActionParameter< bool > preserve_data_format_;
	boost::shared_ptr< ArithmeticFilterAlgo > algo_;
};

bool ActionArithmeticFilterPrivate::validate_parser( Core::ActionContextHandle& context )
{
	// Validate algorithm
	// Create algorithm
	this->algo_ = boost::shared_ptr< ArithmeticFilterAlgo >( new ArithmeticFilterAlgo );

	bool replace = this->replace_.value();

	const std::vector< std::string >& layer_ids = this->layer_ids_.value();
	std::vector< LayerHandle > layers( layer_ids.size() );
	std::string name( "A" );
	for ( size_t i = 0; i < layer_ids.size(); ++i )
	{
		name[ 0 ] = static_cast< char >( 'A' + i );
		this->algo_->find_layer( layer_ids[ i ], layers[ i ] );
		Core::DataBlockHandle input_data_block;
		switch ( layers[ i ]->type() )
		{
		case Core::VolumeType::MASK_E:
			Core::MaskDataBlockManager::Convert( dynamic_cast< MaskLayer* >( layers[ i ].get() )->
				get_mask_volume()->get_mask_data_block(), input_data_block, Core::DataType::UCHAR_E );
			break;
		case Core::VolumeType::DATA_E:
			input_data_block = dynamic_cast< DataLayer* >( layers[ i ].get() )->
				get_data_volume()->get_data_block();
			break;
		}

		std::string error;
		if ( !this->algo_->engine_.add_input_data_block( name, input_data_block, error ) )
		{
			context->report_error( error );
			return false;
		}

		if ( i == 0 && replace )
		{
			this->algo_->lock_for_processing( layers[ i ] );
		}
		else
		{
			this->algo_->lock_for_use( layers[ i ] );
		}
	}

	if ( replace )
	{
		this->algo_->dst_layer_ = layers[ 0 ];
	}
	else if ( this->output_type_.value() == ActionArithmeticFilter::DATA_C )
	{
		this->algo_->create_and_lock_data_layer_from_layer( layers[ 0 ], this->algo_->dst_layer_ );
	}
	else
	{
		this->algo_->create_and_lock_mask_layer_from_layer( layers[ 0 ], this->algo_->dst_layer_ );
	}

	const Core::GridTransform& grid_trans = layers[ 0 ]->get_grid_transform();
	std::string error;
	if ( this->algo_->dst_layer_->type() == Core::VolumeType::MASK_E )
	{
		if( !this->algo_->engine_.add_output_data_block( ActionArithmeticFilter::RESULT_C, 
			grid_trans.get_nx(), grid_trans.get_ny(), grid_trans.get_nz(), 
			Core::DataType::UCHAR_E, error ) )
		{
			context->report_error( error );
		}
	}
	else
	{
		if( !this->algo_->engine_.add_output_data_block( ActionArithmeticFilter::RESULT_C, 
			grid_trans.get_nx(), grid_trans.get_ny(), grid_trans.get_nz(), 
			Core::DataType::FLOAT_E, error ) )
		{
			context->report_error( error );
		}
	}

	this->algo_->engine_.add_expressions( this->expressions_.value() );

	
	if( !this->algo_->engine_.parse_and_validate( error ) )
	{
		context->report_error( error );
		// We need algo to go out of scope so that layers are unlocked
		this->algo_.reset();
		return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////
// Class ActionArithmeticFilter
//////////////////////////////////////////////////////////////////////////

const std::string ActionArithmeticFilter::DATA_C( "data" );
const std::string ActionArithmeticFilter::MASK_C( "mask" );
const std::string ActionArithmeticFilter::RESULT_C( "RESULT" );

ActionArithmeticFilter::ActionArithmeticFilter() :
	private_( new ActionArithmeticFilterPrivate )
{
	// Action arguments
	this->add_argument( this->private_->layer_ids_ );
	this->add_argument( this->private_->expressions_ );
	this->add_argument( this->private_->output_type_ );

	// Action options
	this->add_key( this->private_->replace_ );
	this->add_key( this->private_->preserve_data_format_ );
}

bool ActionArithmeticFilter::validate( Core::ActionContextHandle& context )
{
	const std::vector< std::string >& layer_ids = this->private_->layer_ids_.value();
	if ( layer_ids.size() == 0 )
	{
		context->report_error( "No input layers specified" );
		return false;
	}
	
	// Validate layers
	LayerGroupHandle layer_group;
	for ( size_t i = 0; i < layer_ids.size(); ++i )
	{
		// Check for layer existence
		LayerHandle layer = LayerManager::Instance()->get_layer_by_id( layer_ids[ i ] );
		if ( !layer )
		{
			context->report_error( "Layer '" + layer_ids[ i ] + "' doesn't exist" );
			return false;
		}

		// Make sure that all the layers are in the same group
		if ( !layer_group )
		{
			layer_group = layer->get_layer_group();
		}
		else if ( layer_group != layer->get_layer_group() )
		{
			context->report_error( "Input layers do not belong to the same group" );
			return false;
		}
		
		// Check for layer availability 
		Core::NotifierHandle notifier;
		if ( !LayerManager::CheckLayerAvailability( layer_ids[ i ], 
			i == 0 && this->private_->replace_.value(), notifier ) )
		{
			context->report_need_resource( notifier );
			return false;
		}
	}

	// Validate parser inputs, outputs, and expression
	if( !this->private_->validate_parser( context ) )
	{
		return false;
	}
	
	// Validation successful
	return true;
}

bool ActionArithmeticFilter::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	this->private_->algo_->connect_abort( this->private_->algo_->dst_layer_ );

	// Connect ArrayMathEngine progress signal to output layer progress signal
	this->private_->algo_->engine_.update_progress_signal_.connect(
		boost::bind( &Layer::update_progress, this->private_->algo_->dst_layer_, _1, 0.0, 1.0 ) );

	// Return the ids of the destination layer.
	result = Core::ActionResultHandle( 
		new Core::ActionResult( this->private_->algo_->dst_layer_->get_layer_id() ) );

	// Start the filter.
	Core::Runnable::Start( this->private_->algo_ );

	// We need algo to go out of scope once it finishes so that layers are unlocked
	this->private_->algo_.reset();

	return true;
}

void ActionArithmeticFilter::Dispatch( Core::ActionContextHandle context, 
									  const std::vector< std::string >& layer_ids, 
									  const std::string& expressions, 
									  const std::string& output_type, 
									  bool replace, bool preserve_data_format )
{
	ActionArithmeticFilter* action = new ActionArithmeticFilter;
	action->private_->layer_ids_.set_value( layer_ids );
	action->private_->expressions_.set_value( expressions );
	action->private_->output_type_.set_value( output_type );
	action->private_->replace_.set_value( replace );
	action->private_->preserve_data_format_.set_value( preserve_data_format );

	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

} // end namespace Seg3D
