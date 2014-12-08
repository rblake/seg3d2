/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2014 Scientific Computing and Imaging Institute,
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

#include <Core/Action/ActionFactory.h>

#include <Application/Provenance/Provenance.h>
#include <Application/ProjectManager/ProjectManager.h>
#include <Application/ToolManager/ToolManager.h>
#include <Application/Layer/LayerManager.h>

#include <rbf/Application/Tools/Actions/ActionRBF.h>


CORE_REGISTER_ACTION( Seg3D, RBF )

namespace Seg3D
{

// TODO really worth having a private class here?
class ActionRBFPrivate
{
public:
  std::string target_layer_id_;
  SandboxID sandbox_;
};

ActionRBF::ActionRBF() :
private_( new ActionRBFPrivate )
{
  this->add_parameter( this->private_->sandbox_ );
}

bool ActionRBF::validate( Core::ActionContextHandle& context )
{
  return false;
}

bool ActionRBF::run( Core::ActionContextHandle& context, Core::ActionResultHandle& result )
{
  return false;
}

void ActionRBF::Dispatch( Core::ActionContextHandle context)
{
  ActionRBF* action = new ActionRBF;
  
  Core::ActionDispatcher::PostAction( Core::ActionHandle( action ), context );
}

} // end namespace Seg3D
