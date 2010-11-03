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

// ITK includes
#include <itkConfidenceConnectedImageFilter.h>

// Application includes
#include <Application/LayerManager/LayerManager.h>
#include <Application/StatusBar/StatusBar.h>
#include <Application/Filters/ITKFilter.h>
#include <Application/Filters/Actions/ActionConfidenceConnectedFilter.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
// NOTE: Registration needs to be done outside of any namespace
CORE_REGISTER_ACTION( Seg3D, ConfidenceConnectedFilter )

namespace Seg3D
{

bool ActionConfidenceConnectedFilter::validate( Core::ActionContextHandle& context )
{
	// Check for layer existance and type information
	std::string error;
	if ( ! LayerManager::CheckLayerExistanceAndType( this->target_layer_.value(), 
		Core::VolumeType::DATA_E, error ) )
	{
		context->report_error( error );
		return false;
	}
	
	// Check for layer availability 
	Core::NotifierHandle notifier;
	if ( ! LayerManager::CheckLayerAvailabilityForProcessing( this->target_layer_.value(), 
		notifier ) )
	{
		context->report_need_resource( notifier );
		return false;
	}

	if ( this->seeds_.value().size() == 0 )
	{
		context->report_error( "There needs to be at least one seed point." );
		return false;
	}
		
	if ( this->iterations_.value() < 1 )
	{
		context->report_error( "The number of iterations must be at least 1." );
		return false;
	}
	
	// If the number of iterations is lower than one, we cannot run the filter
	if( this->multiplier_.value() <= 0.0 )
	{
		context->report_error( "The multiplier needs to be larger than zero." );
		return false;
	}
	
	// Validation successful
	return true;
}

// ALGORITHM CLASS
// This class does the actual work and is run on a separate thread.
// NOTE: The separation of the algorithm into a private class is for the purpose of running the
// filter on a separate thread.

class ConfidenceConnectedFilterAlgo : public ITKFilter
{

public:
	LayerHandle src_layer_;
	LayerHandle dst_layer_;

	std::vector< Core::Point > seeds_;
	unsigned int iterations_;
	double multiplier_;
	
public:
	// RUN:
	// Implementation of run of the Runnable base class, this function is called when the thread
	// is launched.
	
	// NOTE: The macro needs a data type to select which version to run. This needs to be
	// a member variable of the algorithm class.
	SCI_BEGIN_ITK_RUN()
	{
		// Define the type of filter that we use.
		typedef itk::ConfidenceConnectedImageFilter< 
			FLOAT_IMAGE_TYPE, FLOAT_IMAGE_TYPE > filter_type;

		// Retrieve the image as an itk image from the underlying data structure
		// NOTE: This only does wrapping and does not regenerate the data.
		Core::ITKImageDataT<float>::Handle input_image; 
		this->get_itk_image_from_layer<float>( this->src_layer_, input_image );
				
		// Create a new ITK filter instantiation.		
		filter_type::Pointer filter = filter_type::New();

		// Relay abort and progress information to the layer that is executing the filter.
		this->observe_itk_filter( filter, this->dst_layer_ );

		// Setup the filter parameters that we do not want to change.
		filter->SetInput( input_image->get_image() );
		filter->SetNumberOfIterations( this->iterations_ );
		filter->SetMultiplier( this->multiplier_ );

		for ( size_t i = 0; i < this->seeds_.size(); ++i )
		{
			filter_type::IndexType index;
			index[ 0 ] = static_cast< int >( this->seeds_[ i ][ 0 ] );
			index[ 1 ] = static_cast< int >( this->seeds_[ i ][ 1 ] );
			index[ 2 ] = static_cast< int >( this->seeds_[ i ][ 2 ] );
			filter->AddSeed( index );
		}
		
		// Run the actual ITK filter.
		// This needs to be in a try/catch statement as certain filters throw exceptions when they
		// are aborted. In that case we will relay a message to the status bar for information.
		try 
		{ 
			filter->Update(); 
		} 
		catch ( ... ) 
		{
			if ( this->check_abort() )
			{
				this->report_error( "Filter was aborted." );
				return;
			}

			this->report_error( "Internal error." );
			return;
		}

		// As ITK filters generate an inconsistent abort behavior, we record our own abort flag
		// This one is set when the abort button is pressed and an abort is sent to ITK.
		if ( this->check_abort() ) return;
		
		this->insert_itk_image_into_layer( this->dst_layer_, filter->GetOutput() );	
	}
	SCI_END_ITK_RUN()
	
	// GET_FITLER_NAME:
	// The name of the filter, this information is used for generating new layer labels.
	virtual std::string get_filter_name() const
	{
		return "ConfidenceConnected Filter";
	}

	// GET_LAYER_PREFIX:
	// This function returns the name of the filter. The latter is prepended to the new layer name, 
	// when a new layer is generated. 
	virtual std::string get_layer_prefix() const
	{
		return "ConfidenceConnected";	
	}
};


bool ActionConfidenceConnectedFilter::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	// Create algorithm
	boost::shared_ptr<ConfidenceConnectedFilterAlgo> algo( new ConfidenceConnectedFilterAlgo );

	// Copy the parameters over to the algorithm that runs the filter
	algo->iterations_ = this->iterations_.value();
	algo->multiplier_ = this->multiplier_.value();

	// Find the handle to the layer
	algo->find_layer( this->target_layer_.value(), algo->src_layer_ );

	// Check the seed points against the source layer dimensions
	Core::GridTransform grid_trans = algo->src_layer_->get_grid_transform();
	Core::Transform inverse_trans = grid_trans.get_inverse();
	int nx = static_cast< int >( grid_trans.get_nx() );
	int ny = static_cast< int >( grid_trans.get_ny() );
	int nz = static_cast< int >( grid_trans.get_nz() );
	const std::vector< Core::Point > seeds = this->seeds_.value();
	for ( size_t i = 0; i < seeds.size(); ++i )
	{
		Core::Point seed = inverse_trans * seeds[ i ];
		seed[ 0 ] = Core::Round( seed[ 0 ] );
		seed[ 1 ] = Core::Round( seed[ 1 ] );
		seed[ 2 ] = Core::Round( seed[ 2 ] );
		if ( seed[ 0 ] >= 0 && seed[ 0 ] < nx &&
			seed[ 1 ] >= 0 && seed[ 1 ] < ny &&
			seed[ 2 ] >= 0 && seed[ 2 ] < nz )
		{
			algo->seeds_.push_back( seed );
		}
	}
	if ( algo->seeds_.size() == 0 )
	{
		context->report_error( "All seed points are out of the volume boundary" );
		return false;
	}
	
	// Lock the src layer, so it cannot be used else where
	algo->lock_for_use( algo->src_layer_ );
	
	// Create the destination layer, which will show progress
	algo->create_and_lock_mask_layer_from_layer( algo->src_layer_, algo->dst_layer_ );

	// Return the id of the destination layer.
	result = Core::ActionResultHandle( new Core::ActionResult( algo->dst_layer_->get_layer_id() ) );

	// Build the undo-redo record
	algo->create_undo_redo_record( context, this->shared_from_this() );

	// Start the filter on a separate thread.
	Core::Runnable::Start( algo );

	return true;
}

void ActionConfidenceConnectedFilter::Dispatch( Core::ActionContextHandle context, 
		std::string target_layer, const std::vector< Core::Point >& seeds,
		unsigned int iterations, double multiplier )
{	
	// Create a new action
	ActionConfidenceConnectedFilter* action = new ActionConfidenceConnectedFilter;

	// Setup the parameters
	action->target_layer_.value() = target_layer;
	action->seeds_.value() = seeds;
	action->iterations_.value() = iterations;
	action->multiplier_.value() = multiplier;

	// Dispatch action to underlying engine
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}
	
} // end namespace Seg3D
