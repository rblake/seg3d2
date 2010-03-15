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

MESSAGE(STATUS "Configuring External Packages")

##################################################
# ADDING ZLIB
##################################################

MESSAGE(STATUS "Configuring zlib")
ADD_SUBDIRECTORY(zlib)

##################################################
# ADDING MATLABIO
##################################################

MESSAGE(STATUS "Configuring MatlabIO")
ADD_SUBDIRECTORY(MatlabIO)

##################################################
# ADDING LIBPNG
##################################################

MESSAGE(STATUS "Configuring libpng")
ADD_SUBDIRECTORY(libpng)

##################################################
# ADDING TEEM
##################################################

# Teem depends on libpng and zlib, so make sure 
# that these projects are included before Teem.

MESSAGE(STATUS "Configuring teem")
ADD_SUBDIRECTORY(teem)

##################################################
# ADDING GLEW
##################################################

MESSAGE(STATUS "Configuring glew")
ADD_SUBDIRECTORY(glew)

##################################################
# ADDING SLIVR
##################################################

# FOR NOW SLIVR IS DISABLED AS IT DOES NOT COMPILE
# WITH GLEWMX - WHICH IS THE CORRECT GLEW SETTING

MESSAGE(STATUS "Configuring slivr")
ADD_SUBDIRECTORY(slivr)

##################################################
# ADDING LIBXML2
##################################################

MESSAGE(STATUS "Configuring libxml2")
ADD_SUBDIRECTORY(libxml2)

##################################################
# ADDING FREETYPE
##################################################

MESSAGE(STATUS "Configuring freetype")
ADD_SUBDIRECTORY(freetype)

##################################################
# ADDING BOOST
##################################################

MESSAGE(STATUS "Configuring boost")
ADD_SUBDIRECTORY(boost)

##################################################
# ADDING ITK
##################################################

MESSAGE(STATUS "Configuring insight toolkit")

# Force ITK not to install any libraries, header
# files etc. Because we don't need that stuff

SET(ITK_INSTALL_NO_DEVELOPMENT 1)
SET(ITK_INSTALL_NO_RUNTIME 1)
SET(ITK_INSTALL_NO_DOCUMENTATION 1)
SET(ITK_INSTALL_NO_LIBRARIES 1)

ADD_SUBDIRECTORY(itk)


