#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
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

OPTION(USE_PRECOMPILED_HEADERS "Use precompiled headers to speed up compilation" OFF)

MACRO(CORE_ADD_LIBRARY name)

  ADD_LIBRARY( ${name} STATIC ${ARGN})

  IF(USE_PRECOMPILED_HEADERS)
    IF(${CMAKE_GENERATOR} MATCHES "Xcode")
	  FILE( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Precompiled )
  
      SET_TARGET_PROPERTIES( ${name} PROPERTIES 
	    XCODE_ATTRIBUTE_SHARED_PRECOMPS_DIR ${CMAKE_BINARY_DIR}/Precompiled
	    XCODE_ATTRIBUTE_GCC_PREFIX_HEADER ${CMAKE_SOURCE_DIR}/Configuration/PrefixHeader.h
	    XCODE_ATTRIBUTE_GCC_PRECOMPILE_PREFIX_HEADER YES
	    XCODE_ATTRIBUTE_PRECOMPS_INCLUDE_HEADERS_FROM_BUILT_PRODUCTS_DIR NO)
    ENDIF(${CMAKE_GENERATOR} MATCHES "Xcode")
  ENDIF(USE_PRECOMPILED_HEADERS)


ENDMACRO(CORE_ADD_LIBRARY)

