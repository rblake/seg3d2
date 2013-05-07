#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2013 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software. 
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#

##################################################
# Set sources
##################################################

SET(APPLICATION_BACKSCATTERRECONSTRUCTION_SRCS
  CalibrationTool.cc
  CalibrationTool.h
  ReconstructionTool.cc
  ReconstructionTool.h
  ViewerTool.cc
  ViewerTool.h
)

SET(APPLICATION_BACKSCATTERRECONSTRUCTION_ACTIONS_SRCS
  Actions/ActionCalibrationView.cc
  Actions/ActionCalibrationView.h
  Actions/ActionReconstructionView.cc
  Actions/ActionReconstructionView.h
  Actions/ActionVisualizationView.cc
  Actions/ActionVisualizationView.h
)

IF(BUILD_WITH_PYTHON)
  GENERATE_ACTION_PYTHON_WRAPPER(PYTHON_WRAPPER
    Application_BackscatterReconstruction ${APPLICATION_BACKSCATTERRECONSTRUCTION_ACTIONS_SRCS})
  SET(APPLICATION_BACKSCATTERRECONSTRUCTION_NOREGISTER_SRCS ${PYTHON_WRAPPER})
ENDIF(BUILD_WITH_PYTHON)

CORE_ADD_LIBRARY(Application_BackscatterReconstruction
  ${APPLICATION_BACKSCATTERRECONSTRUCTION_NOREGISTER_SRCS}
  ${APPLICATION_BACKSCATTERRECONSTRUCTION_SRCS}
  ${APPLICATION_BACKSCATTERRECONSTRUCTION_ACTIONS_SRCS} )
            

TARGET_LINK_LIBRARIES(Application_BackscatterReconstruction
  Core_Utils
  Core_Application
  Core_Interface
  Core_Action
  Core_State
  Core_Geometry
  Application_Clipboard
  Application_Filters
  Application_Tool
  Application_ToolManager
  Application_ViewerManager
  ${SCI_BOOST_LIBRARY}
  ${SCI_ITK_LIBRARIES})


# Register action classes
REGISTER_LIBRARY_AND_CLASSES(Application_BackscatterReconstruction
  ${APPLICATION_BACKSCATTERRECONSTRUCTION_SRCS}
  ${APPLICATION_BACKSCATTERRECONSTRUCTION_ACTIONS_SRCS})

ADD_SUBDIRECTORY(Algorithm)
