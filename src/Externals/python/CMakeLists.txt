cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(Python)

##########################################################################
# Version information                                                    #
##########################################################################
FILE(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/Include/patchlevel.h PY_MAJOR_VERSIONSTR
  REGEX "#define[ \t]+PY_MAJOR_VERSION[ \t]+[0-9]+")
STRING(REGEX MATCH "[0-9]+" PY_MAJOR_VERSIONSTR ${PY_MAJOR_VERSIONSTR})
STRING(LENGTH ${PY_MAJOR_VERSIONSTR} PY_MAJOR_VERSIONSTR_LEN)

FILE(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/Include/patchlevel.h PY_MINOR_VERSIONSTR
  REGEX "#define[ \t]+PY_MINOR_VERSION[ \t]+[0-9]+")
STRING(REGEX MATCH "[0-9]+" PY_MINOR_VERSIONSTR ${PY_MINOR_VERSIONSTR})
STRING(LENGTH ${PY_MINOR_VERSIONSTR} PY_MINOR_VERSIONSTR_LEN)

FILE(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/Include/patchlevel.h PY_MICRO_VERSIONSTR
  REGEX "#define[ \t]+PY_MICRO_VERSION[ \t]+[0-9]+")
STRING(REGEX MATCH "[0-9]+" PY_MICRO_VERSIONSTR ${PY_MICRO_VERSIONSTR})
STRING(LENGTH ${PY_MICRO_VERSIONSTR} PY_MICRO_VERSIONSTR_LEN)

IF(PY_MAJOR_VERSIONSTR_LEN AND PY_MINOR_VERSIONSTR_LEN AND PY_MICRO_VERSIONSTR_LEN)
  MATH(EXPR PY_MAJOR "${PY_MAJOR_VERSIONSTR}")
  MATH(EXPR PY_MINOR "${PY_MINOR_VERSIONSTR}")
  MATH(EXPR PY_PATCH "${PY_MICRO_VERSIONSTR}")
  SET(PYTHON_VERSION "${PY_MAJOR}.${PY_MINOR}.${PY_PATCH}")
  MESSAGE(STATUS "Python version ${PYTHON_VERSION}")
ELSE(PY_MAJOR_VERSIONSTR_LEN AND PY_MINOR_VERSIONSTR_LEN AND PY_MICRO_VERSIONSTR_LEN)
  MESSAGE(FATAL_ERROR "Unable to parse Python version from ${CMAKE_CURRENT_SOURCE_DIR}/Include/patchlevel.h")
ENDIF(PY_MAJOR_VERSIONSTR_LEN AND PY_MINOR_VERSIONSTR_LEN AND PY_MICRO_VERSIONSTR_LEN)

# user configurable options

option(WITH_TSC "profile with the Pentium timestamp counter" 0)
option(ENABLE_IPV6 "Enable IP v6" 1)
option(WITH_DOC_STRINGS "Enable if you want documentation strings in extension modules" 1)
option(Py_USING_UNICODE "Enable unicode support" 1)
option(Py_DEBUG "Build debugging interpreter" 0)
option(WITH_PYMALLOC "Define if you want to compile in Python-specific mallocs" 1)
option(HAVE_LARGEFILE_SUPPORT "Large file support "1)
option(WITH_THREAD "compile in rudimentary thread support" 1)

MARK_AS_ADVANCED(FORCE WITH_TSC ENABLE_IPV6 WITH_DOC_STRINGS Py_USING_UNICODE 
				Py_DEBUG WITH_PYMALLOC HAVE_LARGEFILE_SUPPORT WITH_THREAD)

# do the configure checks

include(ConfigureChecks.cmake)


set(extraLibs "")

set(THREAD_SRC)
if(WITH_THREAD)
  if(HAVE_LIBPTHREAD)
    set(THREAD_SRC Python/thread.c)
    set(_REENTRANT 1)
    set(extraLibs ${extraLibs} ${HAVE_LIBPTHREAD})
  endif(HAVE_LIBPTHREAD)
  if(WIN32)
    set(THREAD_SRC Python/thread.c Python/thread_nt.h)
  endif(WIN32)
endif(WITH_THREAD)

if(HAVE_LIBM)
  set(extraLibs ${extraLibs} ${HAVE_LIBM})
endif(HAVE_LIBM)

if(HAVE_LIBUTIL)
  set(extraLibs ${extraLibs} ${HAVE_LIBUTIL})
endif(HAVE_LIBUTIL)

set(DYNLOAD_SRC)
if(HAVE_DYNAMIC_LOADING)
  if(HAVE_LIBDL)
    set(extraLibs ${extraLibs} ${HAVE_LIBDL})
    set(DYNLOAD_SRC ${DYNLOAD_SRC} Python/dynload_shlib.c)
  endif(HAVE_LIBDL)
  if(WIN32)
    set(DYNLOAD_SRC ${DYNLOAD_SRC} Python/dynload_win.c PC/import_nt.c PC/dl_nt.c)
  endif(WIN32)
#  set(BUILD_LIBS_SHARED_DEFAULT TRUE)
else(HAVE_DYNAMIC_LOADING)
#  set(BUILD_LIBS_SHARED_DEFAULT FALSE)
endif(HAVE_DYNAMIC_LOADING)

set(BUILD_LIBS_SHARED_DEFAULT FALSE)

#option(PYTHON_BUILD_LIB_SHARED "Build libpython shared" ${BUILD_LIBS_SHARED_DEFAULT})
#if(PYTHON_BUILD_LIB_SHARED)
#  set(LIB_BUILD_TYPE SHARED)
#else(PYTHON_BUILD_LIB_SHARED)
  set(LIB_BUILD_TYPE STATIC)
#endif(PYTHON_BUILD_LIB_SHARED)

set(MACHDEP_SRCS )   #TODO: configure

set(SIGNAL_SRCS )    #TODO: configure

set(MOD_SRCS )


#####################################
#
# setup modules
#
#####################################

set(STATIC_MODULES)
set(DYNAMIC_MODULES)

macro(ADD_PYTHON_MODULE _NAME _STATIC_SRCS)
  #option(MODULE_${_NAME}_ENABLE "Add module ${_NAME}" TRUE)
  #MARK_AS_ADVANCED(FORCE MODULE_${_NAME}_ENABLE)
  SET(MODULE_${_NAME}_ENABLE TRUE)

  if(MODULE_${_NAME}_ENABLE)
    #option(MODULE_${_NAME}_SHARED "Add module ${_NAME} shared" ${BUILD_LIBS_SHARED_DEFAULT})
	SET(MODULE_${_NAME}_SHARED FALSE)
    if(MODULE_${_NAME}_SHARED)
      message(STATUS "Adding ${_NAME} as shared module")
      set(DYNAMIC_MODULES ${DYNAMIC_MODULES} ${_NAME})
      add_library(${_NAME} MODULE ${ARGN})
      set_target_properties(${_NAME} PROPERTIES PREFIX "")
    else(MODULE_${_NAME}_SHARED)
      message(STATUS "Adding ${_NAME} as static module")
      set(STATIC_MODULES ${STATIC_MODULES} ${_NAME})
      set(${_STATIC_SRCS} ${${_STATIC_SRCS}} ${ARGN})
      set(PY_STATIC_MODULES_INIT_FUNCTIONS "${PY_STATIC_MODULES_INIT_FUNCTIONS}extern PyObject* PyInit_${_NAME}(void);\n")
      set(PY_STATIC_MODULES_INITTAB        "${PY_STATIC_MODULES_INITTAB}  {\"${_NAME}\", PyInit_${_NAME}},\n")
    endif(MODULE_${_NAME}_SHARED)
  endif(MODULE_${_NAME}_ENABLE)
endmacro(ADD_PYTHON_MODULE)

if(WITH_THREAD)
  add_python_module(_thread  MOD_SRCS Modules/_threadmodule.c)
endif(WITH_THREAD)

if(WIN32)
  set(PY_POSIX_NAME nt)
else(WIN32)
  set(PY_POSIX_NAME posix)
endif(WIN32)

#_io module
add_python_module(_io         	MOD_SRCS Modules/_io/_iomodule.c)
SET(MOD_SRCS ${MOD_SRCS}
			Modules/_io/bufferedio.c
			Modules/_io/bytesio.c
			Modules/_io/fileio.c
			Modules/_io/iobase.c
			Modules/_io/stringio.c
			Modules/_io/textio.c
			Modules/_math.c
			Modules/_time.c)
			
#zlib module
add_python_module(zlib         	MOD_SRCS Modules/zlibmodule.c)
			
#add_python_module(audioop       MOD_SRCS Modules/audioop.c)
			
add_python_module(atexit        MOD_SRCS Modules/atexitmodule.c)
add_python_module(_pickle       MOD_SRCS Modules/_pickle.c)
add_python_module(array         MOD_SRCS Modules/arraymodule.c)
add_python_module(binascii      MOD_SRCS Modules/binascii.c)
add_python_module(cmath         MOD_SRCS Modules/cmathmodule.c)
add_python_module(errno     	MOD_SRCS Modules/errnomodule.c)
add_python_module(math          MOD_SRCS Modules/mathmodule.c)
add_python_module(${PY_POSIX_NAME}  MOD_SRCS Modules/posixmodule.c)
add_python_module(operator      MOD_SRCS Modules/operator.c)
add_python_module(signal    	MOD_SRCS Modules/signalmodule.c)
add_python_module(_md5			MOD_SRCS Modules/md5module.c)
add_python_module(_sha1			MOD_SRCS Modules/sha1module.c)
add_python_module(_sha256		MOD_SRCS Modules/sha256module.c)
add_python_module(_sha512		MOD_SRCS Modules/sha512module.c)
add_python_module(time			MOD_SRCS Modules/timemodule.c)
add_python_module(_locale		MOD_SRCS Modules/_localemodule.c)
add_python_module(_codecs   	MOD_SRCS Modules/_codecsmodule.c)
add_python_module(_weakref  	MOD_SRCS Modules/_weakref.c)
add_python_module(_random       MOD_SRCS Modules/_randommodule.c)
add_python_module(_bisect       MOD_SRCS Modules/_bisectmodule.c)
add_python_module(_heapq        MOD_SRCS Modules/_heapqmodule.c)
add_python_module(_lsprof       MOD_SRCS Modules/_lsprof.c)
add_python_module(itertools     MOD_SRCS Modules/itertoolsmodule.c)
add_python_module(_collections  MOD_SRCS Modules/_collectionsmodule.c)
add_python_module(_symtable 	MOD_SRCS Modules/symtablemodule.c)
add_python_module(mmap 			MOD_SRCS Modules/mmapmodule.c)
add_python_module(_csv          MOD_SRCS Modules/_csv.c)
add_python_module(_sre      	MOD_SRCS Modules/_sre.c)
add_python_module(parser        MOD_SRCS Modules/parsermodule.c)
add_python_module(_struct       MOD_SRCS Modules/_struct.c)
add_python_module(_datetime      MOD_SRCS Modules/_datetimemodule.c)
add_python_module(_functools    MOD_SRCS Modules/_functoolsmodule.c)
add_python_module(_json    		MOD_SRCS Modules/_json.c )
add_python_module(xxsubtype 	MOD_SRCS Modules/xxsubtype.c)
add_python_module(zipimport 	MOD_SRCS Modules/zipimport.c)

IF(WIN32)
  add_python_module(msvcrt		MOD_SRCS PC/msvcrtmodule.c)
  add_python_module(_subprocess	MOD_SRCS PC/_subprocess.c)
  add_python_module(winreg		MOD_SRCS PC/winreg.c)
ELSE(WIN32)
  add_python_module(pwd		MOD_SRCS Modules/pwdmodule.c)
ENDIF(WIN32)

#CJK codecs
add_python_module(_multibytecodec 	MOD_SRCS Modules/cjkcodecs/multibytecodec.c)
add_python_module(_codecs_cn 		MOD_SRCS Modules/cjkcodecs/_codecs_cn.c)
add_python_module(_codecs_hk 		MOD_SRCS Modules/cjkcodecs/_codecs_hk.c)
add_python_module(_codecs_iso2022 	MOD_SRCS Modules/cjkcodecs/_codecs_iso2022.c)
add_python_module(_codecs_jp 		MOD_SRCS Modules/cjkcodecs/_codecs_jp.c)
add_python_module(_codecs_kr 		MOD_SRCS Modules/cjkcodecs/_codecs_kr.c)
add_python_module(_codecs_tw 		MOD_SRCS Modules/cjkcodecs/_codecs_tw.c)

option(BUILD_UNTESTED_MODULES "Build extra modules" FALSE)
MARK_AS_ADVANCED(FORCE BUILD_UNTESTED_MODULES)
if(BUILD_UNTESTED_MODULES)
  add_python_module(_elementtree  MOD_SRCS Modules/_elementtree.c )
  add_python_module(fcntl         MOD_SRCS Modules/fcntlmodule.c )
  add_python_module(grp           MOD_SRCS Modules/grpmodule.c )
endif(BUILD_UNTESTED_MODULES)


file(MAKE_DIRECTORY ${Python_BINARY_DIR}/Modules)
if(WIN32)
  configure_file(PC/pyconfig.h ${Python_BINARY_DIR}/pyconfig.h COPYONLY)
else(WIN32)
  configure_file(pyconfig.h.cmake ${Python_BINARY_DIR}/pyconfig.h @ONLY)
endif(WIN32)
configure_file(Modules/config.c.cmake ${CMAKE_CURRENT_BINARY_DIR}/Modules/modules-config.c @ONLY)
file(WRITE ${Python_BINARY_DIR}/Modules/Setup "#created by cmake, do not edit\n")


#####################################
#
# build python and libpython
# 
#####################################

set(P_SRCS      Parser/acceler.c
				Parser/bitset.c
				Parser/firstsets.c
				Parser/grammar.c
                Parser/grammar1.c
                Parser/listnode.c
				Parser/metagrammar.c
                Parser/node.c
                Parser/parser.c
				Parser/parser.h
                Parser/parsetok.c
                Parser/pgen.c )

set(PARSER_SRCS ${P_SRCS} Parser/myreadline.c Parser/tokenizer.c Parser/tokenizer.h )

set(OBJECT_SRCS Objects/abstract.c
                Objects/boolobject.c
                Objects/bytearrayobject.c
				Objects/bytes_methods.c
				Objects/bytesobject.c
				Objects/capsule.c
                Objects/cellobject.c
                Objects/classobject.c
                Objects/codeobject.c
                Objects/complexobject.c
                Objects/descrobject.c
				Objects/dictobject.c
                Objects/enumobject.c
                Objects/exceptions.c
                Objects/fileobject.c
                Objects/floatobject.c
                Objects/frameobject.c
                Objects/funcobject.c
                Objects/genobject.c
                Objects/iterobject.c
                Objects/listobject.c
                Objects/longobject.c
				Objects/memoryobject.c
                Objects/methodobject.c
                Objects/moduleobject.c
                Objects/object.c
                Objects/obmalloc.c
                Objects/rangeobject.c
                Objects/setobject.c
                Objects/sliceobject.c
                Objects/structseq.c
                Objects/tupleobject.c
                Objects/typeobject.c
                Objects/weakrefobject.c
   )

if(Py_USING_UNICODE)
   set(OBJECT_SRCS ${OBJECT_SRCS} 
                   Objects/unicodeobject.c 
                   Objects/unicodectype.c)
endif(Py_USING_UNICODE)

set(PYTHON_SRCS Python/Python-ast.c 
				Python/_warnings.c
                Python/asdl.c 
                Python/ast.c 
                Python/bltinmodule.c 
                Python/ceval.c 
                Python/compile.c 
                Python/codecs.c 
				Python/dtoa.c
                Python/errors.c 
				Python/fileutils.c
                Python/formatter_unicode.c
                Python/frozen.c 
                Python/frozenmain.c 
                Python/future.c 
                Python/getargs.c 
                Python/getcompiler.c 
                Python/getcopyright.c 
                Python/getopt.c 
                Python/getplatform.c 
                Python/getversion.c 
                Python/graminit.c 
                Python/import.c 
                Python/importdl.c 
				Python/importdl.h
                Python/marshal.c 
                Python/modsupport.c 
                Python/mystrtoul.c 
                Python/mysnprintf.c 
                Python/peephole.c
                Python/pyarena.c 
				Python/pyctype.c
                Python/pyfpe.c 
				Python/pymath.c
                Python/pystate.c 
                Python/pystrcmp.c
                Python/pystrtod.c 
                Python/pythonrun.c
				Python/pytime.c
                Python/structmember.c 
                Python/symtable.c 
                Python/sysmodule.c 
                Python/traceback.c 
                ${MACHDEP_SRCS}
                ${THREAD_SRC}
				${DYNLOAD_SRC}
                )

set(MODULE_SRCS ${CMAKE_CURRENT_BINARY_DIR}/Modules/modules-config.c
				Modules/getbuildinfo.c 
                Modules/main.c 
                Modules/gcmodule.c
				Modules/rotatingtree.c
                ${MOD_SRCS})
				
IF(WIN32)
	SET(MODULE_SRCS ${MODULE_SRCS}
					#PC/config.c
					PC/getpathp.c)
#	set_source_files_properties(PC/getpathp.c PROPERTIES 
#					COMPILE_FLAGS  "-DLANDMARK=L\\\"python\\\\lib\\\\os.py\\\" -DPYTHONPATH=L\\\".\\\\python\\\\lib\\\" ")
ELSE(WIN32)
	SET(MODULE_SRCS ${MODULE_SRCS}
					Modules/getpath.c)
#	set_source_files_properties(Modules/getpath.c PROPERTIES COMPILE_FLAGS  -DVERSION=\\\"${PY_MAJOR}.${PY_MINOR}\\\")
	set_source_files_properties(Python/getplatform.c PROPERTIES COMPILE_FLAGS -DPLATFORM="\\\\\"${PY_PLATFORM}\\\\\"")
ENDIF(WIN32)				


include_directories(Include ${Python_BINARY_DIR} ${SCI_ZLIB_INCLUDE})
if(WIN32)
  include_directories(PC Python)
  add_definitions(-DWIN32 -D_WINDOWS)
endif(WIN32)

add_definitions(-DPy_BUILD_CORE)
ADD_DEFINITIONS(-DPy_NO_ENABLE_SHARED)

IF(SCI_ZLIB_MANGLE)
	ADD_DEFINITIONS(-DZ_PREFIX)
ENDIF(SCI_ZLIB_MANGLE)

#if(NOT HAVE_HYPOT)
#  add_definitions(-DWITHOUT_COMPLEX)
#endif(NOT HAVE_HYPOT)

add_library(python  STATIC
              ${PARSER_SRCS}
              ${OBJECT_SRCS}
              ${PYTHON_SRCS}
              ${MODULE_SRCS}
              ${SIGNAL_SRCS}
              )

target_link_libraries(python ${extraLibs} ${SCI_ZLIB_LIBRARY})
#set_target_properties(pythonLib PROPERTIES OUTPUT_NAME python${PY_MAJOR}.${PY_MINOR})
IF(WIN32)
	SET_TARGET_PROPERTIES(python
						PROPERTIES
						COMPILE_FLAGS "-W3 -wd4267 -wd4244")
	FILE(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/pythonlib.zip" PYTHONLIB_SOURCE_PATH)
	FILE(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/${CMAKE_CFG_INTDIR}/lib/" PYTHONLIB_DST_PATH)
	ADD_CUSTOM_COMMAND(TARGET python 
		   POST_BUILD
		   COMMAND if not exist "${PYTHONLIB_DST_PATH}" md "${PYTHONLIB_DST_PATH}"
		   COMMAND copy /Y "${PYTHONLIB_SOURCE_PATH}" "${PYTHONLIB_DST_PATH}")
ELSE(WIN32)
	SET(PYTHONLIB_DST_PATH "${CMAKE_BINARY_DIR}/${CMAKE_CFG_INTDIR}/lib/")
	ADD_CUSTOM_COMMAND(TARGET python 
		   POST_BUILD
		   COMMAND [ -d ${PYTHONLIB_DST_PATH} ] || mkdir "${PYTHONLIB_DST_PATH}" 
		   COMMAND cp -rf "${CMAKE_CURRENT_SOURCE_DIR}/pythonlib.zip" "${PYTHONLIB_DST_PATH}")		   
ENDIF(WIN32)

SET(PYTHON_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "" FORCE)
SET(SCI_PYTHON_LIBRARY python CACHE INTERNAL "Python library." FORCE)
SET(SCI_PYTHON_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/Include ${Python_BINARY_DIR} CACHE INTERNAL "Python include directory." FORCE)
