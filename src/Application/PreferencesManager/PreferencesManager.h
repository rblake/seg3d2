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

#ifndef APPLICATION_PREFERENCESMANAGER_PREFERENCESMANAGER_H
#define APPLICATION_PREFERENCESMANAGER_PREFERENCESMANAGER_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif 



// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

// Core includes
#include <Core/Utils/StringUtil.h>
#include <Core/Utils/Singleton.h>
#include <Core/Utils/Log.h>
#include <Core/Utils/Exception.h>
#include <Core/State/StateHandler.h>


namespace Seg3D
{

// Forward declaration
class PreferencesManager;

// Class definition
class PreferencesManager : public Core::StateHandler
{
	CORE_SINGLETON( PreferencesManager );

	// -- Constructor/Destructor --
private:
	PreferencesManager();
	virtual ~PreferencesManager();
	
	// -- State variables --
public:
	//General Preferences
	Core::StateStringHandle project_path_state_;
	Core::StateBoolHandle considate_project_state_;
	Core::StateBoolHandle full_screen_on_startup_state_;
	
	//Viewers Preferences
	Core::StateOptionHandle default_viewer_mode_state_;
	Core::StateIntHandle grid_size_state_;
	Core::StateOptionHandle background_color_state_;
	Core::StateBoolHandle show_slice_number_state_;
	Core::StateOptionHandle naming_convention_state_;
	
	//Layers Preferences
	Core::StateRangedDoubleHandle default_layer_opacity_state_;
	Core::StateOptionHandle default_mask_fill_state_;
	Core::StateOptionHandle default_mask_border_state_;
	std::vector< Core::StateColorHandle > color_states_;
	Core::StateColorHandle color_01_state_;
	Core::StateColorHandle color_02_state_;
	Core::StateColorHandle color_03_state_;
	Core::StateColorHandle color_04_state_;
	Core::StateColorHandle color_05_state_;
	Core::StateColorHandle color_06_state_;
	Core::StateColorHandle color_07_state_;
	Core::StateColorHandle color_08_state_;
	Core::StateColorHandle color_09_state_;
	Core::StateColorHandle color_10_state_;
	Core::StateColorHandle color_11_state_;
	Core::StateColorHandle color_12_state_;

	//Interface Controls Preferences
		
	//Sidebars Preferences
	Core::StateBoolHandle show_tools_bar_state_;
	Core::StateBoolHandle show_layermanager_bar_state_;
	Core::StateBoolHandle show_projectmanager_bar_state_;
	Core::StateBoolHandle show_measurement_bar_state_;
	Core::StateBoolHandle show_history_bar_state_;
	
public:
	const std::vector< Core::Color >& get_default_colors()
	{
		return default_colors_;
	}
	
private:
	void initialize_states();
	bool initialize_default_colors();
	
private:
	std::vector< Core::Color > default_colors_;
	
		
};

} // end namespace seg3D

#endif

