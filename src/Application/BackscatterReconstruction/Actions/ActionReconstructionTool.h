/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#ifndef APPLICATION_BACKSCATTERRECONSTRUCTION_ACTIONS_ACTIONRECONSTRUCTIONTOOL_H
#define APPLICATION_BACKSCATTERRECONSTRUCTION_ACTIONS_ACTIONRECONSTRUCTIONTOOL_H


// Core includes
#include <Core/Action/Action.h>
#include <Core/Interface/Interface.h>

#include <Application/Layer/Layer.h>
#include <Application/Layer/DataLayer.h> 
#include <Application/Layer/LayerManager.h>
#include <Application/Layer/LayerAction.h>

#include <boost/smart_ptr.hpp> 

namespace Seg3D
{
  
class ActionReconstructionTool : public LayerAction
{
  
CORE_ACTION(
  CORE_ACTION_TYPE( "ActionReconstructionTool", "" )
  CORE_ACTION_ARGUMENT( "layerid", "The ID of the target layer." )
  CORE_ACTION_OPTIONAL_ARGUMENT( "initialGuessSet", "<none>", "Calibration set." )
  CORE_ACTION_OPTIONAL_ARGUMENT( "outputDir", "", "Algorithm output directory." )
  CORE_ACTION_OPTIONAL_ARGUMENT( "iterations", "3", "Number of iterations to perform." )
  CORE_ACTION_OPTIONAL_ARGUMENT( "xyVoxelSizeScale", "0.5", "Voxel scale in the XY plane." )
  CORE_ACTION_OPTIONAL_ARGUMENT( "zVoxelSizeScale", "0.5", "Voxel scale in the Z plane." )
  CORE_ACTION_OPTIONAL_ARGUMENT( "sandbox", "-1", "The sandbox in which to run the action." )
  CORE_ACTION_ARGUMENT_IS_NONPERSISTENT( "sandbox" )
  CORE_ACTION_CHANGES_PROJECT_DATA()
)
  
  // -- Constructor/Destructor --
public:
  ActionReconstructionTool();
  
  // -- Functions that describe action --
public:
  virtual bool validate( Core::ActionContextHandle& context );
  virtual bool run( Core::ActionContextHandle& context, Core::ActionResultHandle& result );
  virtual void clear_cache();
  
  typedef boost::function< void ( double, double, double ) > callback_type;
  
  // -- Dispatch this action from the interface --
public:
  /// DISPATCH:
  /// Dispatch an action
  static void Dispatch( Core::ActionContextHandle context,
                        std::string target_layer,
                        const std::vector< std::string >& initialGuessSet,
                        std::string outputDir,
                        int iterations,
                        double xyVoxelSizeScale,
                        double zVoxelSizeScale,
                        callback_type callback = 0 );

  static void Abort();
  
private:
  std::string target_layer_;
  std::vector< std::string > initalGuessSet_;
  int iterations_;
  double xyVoxelSizeScale_;
  double zVoxelSizeScale_;
  std::string outputDir_;
  SandboxID sandbox_;
  callback_type callback_;
};
  
} // end namespace Seg3D

#endif
