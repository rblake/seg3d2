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

#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/LayerManager/Actions/ActionLayer.h>

namespace Seg3D
{

bool ActionLayer::cache_layer_handle( Core::ActionContextHandle& context,
	Core::ActionParameter< std::string >& layer_id, 
	Core::ActionCachedHandle< LayerHandle >& layer )
{
	if ( ! layer.handle() )
	{
		layer.handle() = LayerManager::Instance()->get_layer_by_id( layer_id.value() );
		
		if ( ! layer.handle() )
		{
			context->report_error( std::string( "LayerID: '" ) + layer_id.value() + "' is invalid" );
			return false;
		}
	}
	
	return true;
}

bool ActionLayer::cache_group_handle( Core::ActionContextHandle& context,
	Core::ActionParameter< std::string >& group_id, 
	Core::ActionCachedHandle< LayerGroupHandle >& group )
{
	if ( ! group.handle() )
	{
		group.handle() = LayerManager::Instance()->get_layer_group( group_id.value() );
		
		if ( ! group.handle() )
		{
			context->report_error( std::string( "LayerGroupID: '" ) + 
				group_id.value() + "' is invalid" );
			return false;
		}
	}
	
	return true;
}

bool ActionLayer::check_availability( Core::ActionContextHandle& context,
	Core::ActionCachedHandle< LayerHandle >& layer )
{
	Core::ResourceLockHandle resource_lock = layer.handle()->get_resource_lock();
	if ( resource_lock->is_locked() )
	{
		context->report_need_resource( resource_lock );
		return false;
	}
	
	return true;
}

} // end namespace Seg3D
