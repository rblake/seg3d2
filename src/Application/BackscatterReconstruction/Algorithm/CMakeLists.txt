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

#
# library for gui
#
IF(CUDA_FOUND)
  CUDA_ADD_LIBRARY(Application_BackscatterReconstruction_Algorithm # no main.cpp!
    geometry.cpp
    geometry.h
    sart.cpp
    sart.h
    markov.cpp
    markov.cu
    markov.h
    calibration.cpp
    calibration.h
    material.cpp
    material.h
    lmmin.cpp
    lmmin.h
    vec3.h 
    gamma.h
    deviates.h
    ran.h
    nr3.h 
    parallel.h
    gradient.h
    cone.h
    recon_api.cpp
    recon_api.h
    recon_params.h STATIC)

  SET_TARGET_PROPERTIES(Application_BackscatterReconstruction_Algorithm
    PROPERTIES
      COMPILE_DEFINITIONS USE_CUDA=1
  )
ELSE(CUDA_FOUND)
  ADD_LIBRARY(Application_BackscatterReconstruction_Algorithm # no main.cpp!
    STATIC
    geometry.cpp
    geometry.h
    sart.cpp
    sart.h
    markov.cpp
    markov.h
    calibration.cpp
    calibration.h
    material.cpp
    material.h
    lmmin.cpp
    lmmin.h
    vec3.h 
    gamma.h
    deviates.h
    ran.h
    nr3.h 
    parallel.h
    gradient.h
    cone.h
    recon_api.cpp
    recon_api.h
    recon_params.h)
ENDIF(CUDA_FOUND)

TARGET_LINK_LIBRARIES(Application_BackscatterReconstruction_Algorithm
  ${SCI_ITK_LIBRARIES}
)

