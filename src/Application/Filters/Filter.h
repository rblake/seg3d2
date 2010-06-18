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

#ifndef APPLICATION_FILTER_FILTER_H
#define APPLICATION_FILTER_FILTER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

// STL includes
#include <vector>

// Boost includes 
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>

// Core includes
#include <Core/Utils/IntrusiveBase.h>
#include <Core/Utils/Lockable.h>

// Application includes
#include <Application/Layer/LayerFWD.h>

namespace Seg3D
{

// Forward Declaration
class Filter;
typedef boost::intrusive_ptr<Filter> FilterHandle;

// Class definition
class Filter : public Core::IntrusiveBase
{

	// -- Constructor/destructor --
public:

	Filter();
	virtual ~Filter();

public:
	// START:
	// Start the filter on a separate thread
	bool start( LayerHandle layer );
			
protected:

	// RUN_FILTER:
	// The actual implementation of the filter
	virtual bool run_filter( LayerHandle layer ) = 0;

private:
	// Function for running the separate thread
	static void RunFilter( FilterHandle filter, LayerHandle layer );

};

} // end namespace Seg3D

#endif
