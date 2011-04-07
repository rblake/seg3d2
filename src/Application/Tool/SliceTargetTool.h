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

#ifndef APPLICATION_TOOL_SLICETARGETTOOL_H
#define APPLICATION_TOOL_SLICETARGETTOOL_H

// Application includes
#include <Application/Tool/SingleTargetTool.h>

namespace Seg3D
{

// Forward declaration
class SliceTargetToolPrivate;
typedef boost::shared_ptr< SliceTargetToolPrivate > SliceTargetToolPrivateHandle;

// Class definition
class SliceTargetTool : public SingleTargetTool
{

	// -- constructor/destructor --
public:
	SliceTargetTool( int target_volume_type, const std::string& tool_type );	
	virtual ~SliceTargetTool();

	// -- state --
public:
	Core::StateLabeledOptionHandle slice_type_state_;
	Core::StateBoolHandle use_active_viewer_state_;

private:
	SliceTargetToolPrivateHandle private_;

public:
	const static std::string AXIAL_C;
	const static std::string SAGITTAL_C;
	const static std::string CORONAL_C;
};

} // end namespace Seg3D

#endif
