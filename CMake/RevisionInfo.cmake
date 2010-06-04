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

MACRO( GENERATE_REVISION_INFO )

	FIND_PACKAGE( Subversion )
	Subversion_WC_INFO( ${PROJECT_SOURCE_DIR} "SVNINFO" )

	IF( Subversion_FOUND )
		SET( INFO_STRING "#define ${PROJECT_NAME}_REVISION \"${SVNINFO_WC_REVISION}\"\n" )
		SET( INFO_STRING "${INFO_STRING}#define ${PROJECT_NAME}_DATE \"${SVNINFO_WC_LAST_CHANGED_DATE}\"\n" )
		SET( INFO_STRING "${INFO_STRING}#define ${PROJECT_NAME}_REVISIONINFO \"${PROJECT_NAME} Revision: ${SVNINFO_WC_REVISION} Date: ${SVNINFO_WC_LAST_CHANGED_DATE}\"\n\n" )

	ELSE( Subversion_FOUND )
		SET( INFO_STRING "#define ${PROJECT_NAME}_REVISION \"unknown\"\n" )
		SET( INFO_STRING "${INFO_STRING}#define ${PROJECT_NAME}_DATE \"unknown\"\n" )
		SET( INFO_STRING "${INFO_STRING}#define ${PROJECT_NAME}_REVISIONINFO \"${PROJECT_NAME} Revision:unknown Date: unknown\"\n\n" )
	ENDIF( Subversion_FOUND )

	FILE( MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_RevisionInfo" )
	SET( FILENAME "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_RevisionInfo/${PROJECT_NAME}_RevisionInfo.h" )
	FILE( WRITE ${FILENAME} ${INFO_STRING})
	INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_RevisionInfo )

ENDMACRO( GENERATE_REVISION_INFO )
