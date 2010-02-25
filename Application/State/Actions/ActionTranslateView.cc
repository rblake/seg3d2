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

#include <Application/State/Actions/ActionTranslateView.h>
#include <Application/State/StateView2D.h>
#include <Application/State/StateView3D.h>

namespace Seg3D
{

SCI_REGISTER_ACTION(TranslateView)
;

ActionTranslateView::ActionTranslateView()
{
	add_argument( this->stateid_ );
	add_argument( this->offset_ );
}

bool ActionTranslateView::validate( ActionContextHandle &context )
{
	StateBaseHandle state = this->state_weak_handle_.lock();
	if ( !state )
	{
		if ( !( StateEngine::Instance()->get_state( stateid_.value(), state ) ) )
		{
			context->report_error( std::string( "Unknown state variable '" ) + stateid_.value()
			    + "'" );
			return false;
		}

		if ( typeid(*state) != typeid(StateView2D) && typeid(*state) != typeid(StateView3D) )
			{
				context->report_error( std::string( "State variable '" ) + stateid_.value()
				    + "' doesn't support ActionScaleView3D" );
				return false;
			}

			this->state_weak_handle_ = boost::dynamic_pointer_cast< StateViewBase >( state );
		}

		return true;
	}

	bool ActionTranslateView::run( ActionContextHandle& context, ActionResultHandle& result )
	{
		StateViewBaseHandle state = this->state_weak_handle_.lock();

		if ( state )
		{
			state->translate( this->offset_.value() );
			return true;
		}

		return false;
	}

	} // end namespace Seg3D
