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

#ifndef APPLICATION_TOOLS_BOOLEANFILTER_H
#define APPLICATION_TOOLS_BOOLEANFILTER_H

#include <Application/Tool/Tool.h>

namespace Seg3D
{

class BooleanFilter : public Tool
{

SEG3D_TOOL(
SEG3D_TOOL_NAME( "BooleanFilter", "Boolean Algebra with masks" )
SEG3D_TOOL_MENULABEL( "Boolean" )
SEG3D_TOOL_MENU( "filter_mask_to_mask" )
SEG3D_TOOL_SHORTCUT_KEY( "Alt+B" )
SEG3D_TOOL_URL( "http://seg3d.org/" )
)

public:
	BooleanFilter( const std::string& toolid );
	virtual ~BooleanFilter();

	// -- constraint parameters --

	// Constrain viewer to right painting tool when layer is selected
	void target_constraint( std::string layerid );

	// -- activate/deactivate tool --

	virtual void activate();
	virtual void deactivate();
	
private:
	// -- handle updates from layermanager --
	void handle_layers_changed();

	// -- state --
public:
	// Layerid of the target layer
	Core::StateStringHandle mask_a_state_;

	Core::StateStringHandle mask_b_state_;

	Core::StateStringHandle mask_c_state_;

	Core::StateStringHandle mask_d_state_;

	Core::StateOptionHandle example_expressions_state_;

	Core::StateBoolHandle replace_state_;

};

} // end namespace

#endif