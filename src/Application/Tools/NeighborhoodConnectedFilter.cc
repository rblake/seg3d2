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

#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/NeighborhoodConnectedFilter.h>


namespace Seg3D {
  
  // Register the tool into the tool factory
  SCI_REGISTER_TOOL(NeighborhoodConnectedFilter)
  
  
  NeighborhoodConnectedFilter::NeighborhoodConnectedFilter(const std::string& toolid) :
  Tool(toolid)
  {
    // Need to set ranges and default values for all parameters
    add_state("target",target_layer_state_,"<none>","<none>");
    
    // Add constaints, so that when the state changes the right ranges of 
    // parameters are selected
    target_layer_state_->value_changed_signal_.connect(
      boost::bind(&NeighborhoodConnectedFilter::target_constraint,this,_1));
    
    
  }
  
  void
  NeighborhoodConnectedFilter::target_constraint(std::string layerid)
  {
  }
  
  NeighborhoodConnectedFilter::~NeighborhoodConnectedFilter()
  {
    disconnect_all();
  }
  
  void
  NeighborhoodConnectedFilter::activate()
  {
  }
  
  void
  NeighborhoodConnectedFilter::deactivate()
  {
  }
  
} // end namespace Seg3D


