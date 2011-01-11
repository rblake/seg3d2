# Define the srcs for the Information Object Definition
# IOD
SET(IOD_SRCS
  gdcmModule.cxx
  gdcmMacro.cxx
  gdcmModules.cxx
  gdcmMacros.cxx
  gdcmNestedModuleEntries.cxx
  gdcmIODEntry.cxx
  gdcmTableReader.cxx
  gdcmTable.cxx
  gdcmSeries.cxx
  gdcmDefs.cxx
  gdcmDefinedTerms.cxx
  gdcmEnumeratedValues.cxx
  gdcmStudy.cxx
  gdcmPatient.cxx
  gdcmType.cxx
  gdcmUsage.cxx
  gdcmIOD.cxx
  #gdcmXMLDictReader.cxx
  #gdcmXMLPrivateDictReader.cxx
  )

# expat stuff
IF(NOT BUILD_SHARED_LIBS)
  SET_SOURCE_FILES_PROPERTIES(gdcmTableReader.cxx
                              PROPERTIES 
                              COMPILE_FLAGS -DXML_STATIC
                              )
ENDIF(NOT BUILD_SHARED_LIBS)

# Add the include paths
INCLUDE_DIRECTORIES(
  "${GDCM_SOURCE_DIR}/Source/Common"
  "${GDCM_BINARY_DIR}/Source/Common"
  "${GDCM_SOURCE_DIR}/Source/DataStructureAndEncodingDefinition"
  "${GDCM_SOURCE_DIR}/Source/DataDictionary"
  "${GDCM_SOURCE_DIR}/Utilities"
  )

ADD_LIBRARY(gdcmIOD ${IOD_SRCS})
TARGET_LINK_LIBRARIES(gdcmIOD gdcmDSED gdcmCommon ${GDCM_EXPAT_LIBRARIES})
SET_TARGET_PROPERTIES(gdcmIOD PROPERTIES ${GDCM_LIBRARY_PROPERTIES} LINK_INTERFACE_LIBRARIES "gdcmDSED;gdcmCommon")
IF(WIN32)
  SET_TARGET_PROPERTIES(gdcmIOD
                        PROPERTIES
                        COMPILE_FLAGS "-W0 -wd4244 -wd4267 -wd4244"
						LINK_FLAGS "/ignore:4221")
ENDIF(WIN32)

IF(CMAKE_COMPILER_IS_GNUCC)
  SET_TARGET_PROPERTIES(gdcmIOD
                        PROPERTIES
                        COMPILE_FLAGS "-w")
ENDIF(CMAKE_COMPILER_IS_GNUCC)


#-----------------------------------------------------------------------------
FILE(GLOB GDCM_XML_TABLES_GLOB
  "${CMAKE_CURRENT_SOURCE_DIR}/*.xml"
  )

SET(GDCM_XML_TABLES)
FOREACH(filename ${GDCM_XML_TABLES_GLOB})
  SET(GDCM_XML_TABLES "${GDCM_XML_TABLES}\n\"${filename}\",")
ENDFOREACH(filename)

# Populate GDCM_DATA_IMAGES:
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/gdcmTables.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/gdcmTables.h"
  )

