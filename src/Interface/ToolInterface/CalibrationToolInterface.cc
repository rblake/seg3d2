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


// Core includes
#include <Core/Interface/Interface.h>
#include <Core/Utils/Log.h>

// Qt includes
#include <QtGui/QClipboard>
#include <QtGui/QColorDialog>

// Qt Gui Includes
#include <QtUtils/Bridge/QtBridge.h>

// Interface Includes
#include <Interface/ToolInterface/CalibrationToolInterface.h>
#include "ui_CalibrationToolInterface.h"

// Application Includes
#include <Application/Layer/LayerManager.h>
//#include <Application/ViewerManager/Actions/ActionPickPoint.h>

// Core Includes
#include <Core/State/Actions/ActionSetAt.h>

SCI_REGISTER_TOOLINTERFACE( Seg3D, CalibrationToolInterface )

namespace Seg3D
{

class CalibrationToolInterfacePrivate
{
public:
	Ui::CalibrationToolInterface ui_;
	CalibrationToolInterface* interface_;  
};

// constructor
CalibrationToolInterface::CalibrationToolInterface() :
private_( new CalibrationToolInterfacePrivate )
{
  this->private_->interface_ = this;
}

// destructor
CalibrationToolInterface::~CalibrationToolInterface()
{
  this->disconnect_all();
}

bool CalibrationToolInterface::build_widget( QFrame* frame )
{
  return true;
}

} // end namespace Seg3D
