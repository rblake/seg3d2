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

#ifndef INTERFACE_BACKSCATTERRECONSTRUCTION_RECONSTRUCTIONTOOLINTERFACE_H
#define INTERFACE_BACKSCATTERRECONSTRUCTION_RECONSTRUCTIONTOOLINTERFACE_H

// Qt includes
#include <QtCore/QPointer>
#include <QtGui/QTableWidget>

// Boost includes
#include <boost/shared_ptr.hpp>

// Base class of the tool widget include
#include <Interface/Application/ToolWidget.h>

namespace Seg3D
{
  
class ReconstructionToolInterfacePrivate;
typedef boost::shared_ptr< ReconstructionToolInterfacePrivate > ReconstructionToolInterfacePrivateHandle;


class ReconstructionToolInterface : public ToolWidget
{
Q_OBJECT
  
public:
  ReconstructionToolInterface();
  virtual ~ReconstructionToolInterface();
  virtual bool build_widget( QFrame* frame );
	void update_progress_bar( double progress );
	void reset_progress_bar();

private:
  ReconstructionToolInterfacePrivateHandle private_;

public:
  typedef QPointer< ReconstructionToolInterface > qpointer_type;
  
  static void UpdateProgress( qpointer_type qpointer, double progress );
  static void ResetProgress( qpointer_type qpointer );

private Q_SLOTS:
  void triggerSetOutputDir();
  void triggerDataImport();
  void triggerLabelImport();
};
  
} // end namespace Seg3D

#endif
