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
#include <itkOtsuMultipleThresholdsImageFilter.h>

// Application includes
#include <Application/LayerManager/LayerManager.h>
#include <Application/StatusBar/StatusBar.h>
#include <Application/Filters/ITKFilter.h>
#include <Application/Filters/Actions/ActionOtsuThresholdFilter.h>

// REGISTER ACTION:
// Define a function that registers the action. The action also needs to be
// registered in the CMake file.
// NOTE: Registration needs to be done outside of any namespace
CORE_REGISTER_ACTION( Seg3D, OtsuThresholdFilter )

namespace Seg3D
{

bool ActionOtsuThresholdFilter::validate( Core::ActionContextHandle& context )
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
		
	// If the number of iterations is lower than one, we cannot run the filter
	if( this->amount_.value() < 1 )
	{
		context->report_error( "The number of thresholds needs to be at least one." );
		return false;
	}
	
	// Validation successful
	return true;
}

// ALGORITHM CLASS
// This class does the actual work and is run on a separate thread.
// NOTE: The separation of the algorithm into a private class is for the purpose of running the
// filter on a separate thread.

class OtsuThresholdFilterAlgo : public ITKFilter
{

public:
	LayerHandle src_layer_;
	std::vector<LayerHandle> dst_layer_;

	double amount_;
	
public:
	// RUN:
	// Implementation of run of the Runnable base class, this function is called when the thread
	// is launched.
	
	// NOTE: The macro needs a data type to select which version to run. This needs to be
	// a member variable of the algorithm class.
	SCI_BEGIN_TYPED_ITK_RUN( this->src_layer_->get_data_type() )
	{
		// Define the type of filter that we use.
		typedef itk::OtsuMultipleThresholdsImageFilter< 
			TYPED_IMAGE_TYPE, UCHAR_IMAGE_TYPE > filter_type;

		// Retrieve the image as an itk image from the underlying data structure
		// NOTE: This only does wrapping and does not regenerate the data.
		typename Core::ITKImageDataT<VALUE_TYPE>::Handle input_image; 
		this->get_itk_image_from_layer<VALUE_TYPE>( this->src_layer_, input_image );
				
		// Create a new ITK filter instantiation.		
		typename filter_type::Pointer filter = filter_type::New();

		// Relay abort and progress information to the layer that is executing the filter.

		for ( size_t j = 0; j < static_cast<size_t>( this->amount_ + 1 );  j++ )
		{
			this->observe_itk_filter( filter, this->dst_layer_[ j ] );
		}
		
		// Setup the filter parameters that we do not want to change.
		filter->SetInput( input_image->get_image() );
	    filter->SetLabelOffset( 0 );
		filter->SetNumberOfThresholds( this->amount_ );

		// Run the actual ITK filter.
		// This needs to be in a try/catch statement as certain filters throw exceptions when they
		// are aborted. In that case we will relay a message to the status bar for information.
		try 
		{ 
			filter->Update(); 
		} 
		catch ( ... ) 
		{
			this->report_error( "Could not allocate enough memory." );
			return;
		}

		// As ITK filters generate an inconsistent abort behavior, we record our own abort flag
		// This one is set when the abort button is pressed and an abort is sent to ITK.
		if ( this->check_abort() ) return;
		
		for ( size_t j = 0; j < static_cast<size_t>( this->amount_ + 1 );  j++ )
		{		
			this->insert_itk_label_into_mask_layer( this->dst_layer_[ j ], filter->GetOutput(), 
				static_cast<unsigned char>( j ) );	
		}
	}
	SCI_END_TYPED_ITK_RUN()
	
	// GET_FITLER_NAME:
	// The name of the filter, this information is used for generating new layer labels.
	virtual std::string get_filter_name() const
	{
		return "OtsuThreshold Filter";
	}

	// GET_LAYER_PREFIX:
	// This function returns the name of the filter. The latter is prepended to the new layer name, 
	// when a new layer is generated. 
	virtual std::string get_layer_prefix() const
	{
		return "OtsuThreshold";	
	}
};


bool ActionOtsuThresholdFilter::run( Core::ActionContextHandle& context, 
	Core::ActionResultHandle& result )
{
	// Create algorithm
	boost::shared_ptr<OtsuThresholdFilterAlgo> algo( new OtsuThresholdFilterAlgo );

	// Copy the parameters over to the algorithm that runs the filter
	algo->amount_ = this->amount_.value();

	// Find the handle to the layer
	if ( !( algo->find_layer( this->target_layer_.value(), algo->src_layer_ ) ) )
	{
		return false;
	}

	algo->dst_layer_.resize( algo->amount_ + 1 );

	// Lock the src layer, so it cannot be used else where
	algo->lock_for_use( algo->src_layer_ );
	
	// Create the destination layer, which will show progress
	for ( size_t j = 0; j < static_cast<size_t>( algo->amount_ + 1 );  j++ )
	{
		algo->create_and_lock_mask_layer_from_layer( algo->src_layer_, algo->dst_layer_[ j ] );
	}
	
	// Return the id of the destination layer.
	if ( algo->dst_layer_.size() )
	{
		result = Core::ActionResultHandle( new Core::ActionResult( algo->dst_layer_[ 0 ]->get_layer_id() ) );
	}
	
	// Start the filter on a separate thread.
	Core::Runnable::Start( algo );

	return true;
}

void ActionOtsuThresholdFilter::Dispatch( Core::ActionContextHandle context, 
	std::string target_layer, int amount )
{	
	// Create a new action
	ActionOtsuThresholdFilter* action = new ActionOtsuThresholdFilter;

	// Setup the parameters
	action->target_layer_.value() = target_layer;
	action->amount_.value() = amount;

	// Dispatch action to underlying engine
	Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}
	
} // end namespace Seg3D
