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

#ifndef APPLICATION_TOOLS_ACTIONS_ACTIOINFLOODFILL_H
#define APPLICATION_TOOLS_ACTIONS_ACTIOINFLOODFILL_H

#include <Application/LayerManager/Actions/ActionLayer.h>

namespace Seg3D
{

class ActionFloodFillPrivate;
typedef boost::shared_ptr< ActionFloodFillPrivate > ActionFloodFillPrivateHandle;

class FloodFillInfo
{
public:
	std::string target_layer_id_;
	int slice_type_;
	size_t slice_number_;
	std::vector< Core::Point > seeds_;
	std::string data_constraint_layer_id_;
	double min_val_;
	double max_val_;
	bool negative_data_constraint_;
	std::string mask_constraint1_layer_id_;
	bool negative_mask_constraint1_;
	std::string mask_constraint2_layer_id_;
	bool negative_mask_constraint2_;
	bool erase_;
};

class ActionFloodFill : public ActionLayer
{

CORE_ACTION
( 
	CORE_ACTION_TYPE( "Floodfill", "Flood fill the content of a mask slice "
		"starting from seed points." )
	CORE_ACTION_ARGUMENT( "target", "The ID of the target mask layer." )
	CORE_ACTION_ARGUMENT( "slice_type", "The slicing direction." )
	CORE_ACTION_ARGUMENT( "slice_number", "The slice number to be filled." )
	CORE_ACTION_ARGUMENT( "seed_points", "The world coordinates of seed points." )
	CORE_ACTION_KEY( "data_constraint", "<none>", "The ID of data constraint layer." )
	CORE_ACTION_KEY( "min_value", "0", "The minimum data constraint value." )
	CORE_ACTION_KEY( "max_value", "0", "The maximum data constraint value." )
	CORE_ACTION_KEY( "negative_data_constraint", "false", "Whether to negate the data constraint." )
	CORE_ACTION_KEY( "mask_constraint1", "<none>", "The ID of first mask constraint layer." )
	CORE_ACTION_KEY( "negative_mask_constraint1", "false", "Whether to negate the first mask constraint." )
	CORE_ACTION_KEY( "mask_constraint2", "<none>", "The ID of second mask constraint layer." )
	CORE_ACTION_KEY( "negative_mask_constraint2", "false", "Whether to negate the second mask constraint." )
	CORE_ACTION_KEY( "erase", "false", "Whether to erase instead of fill." )
)

public:
	ActionFloodFill();
	virtual ~ActionFloodFill();

	// VALIDATE:
	// Validate the action
	virtual bool validate( Core::ActionContextHandle& context );

	// RUN:
	// Run the action
	virtual bool run( Core::ActionContextHandle& context, Core::ActionResultHandle& result );

	// CLEAR_CACHE:
	// Clear any intermediate objects that were created between dispatch and run.
	virtual void clear_cache();

private:
	ActionFloodFillPrivateHandle private_;

public:
	// DISPATCH:
	// Dispatch the action.
	static void Dispatch( Core::ActionContextHandle context, const FloodFillInfo& params );
};

} // end namespace Seg3D

#endif