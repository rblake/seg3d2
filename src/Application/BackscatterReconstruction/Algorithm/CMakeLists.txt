#cmake_minimum_required(VERSION 2.8)
#cmake_policy(SET CMP0015 NEW)

#PROJECT(backscatter)

# CUDA
#FIND_PACKAGE(CUDA REQUIRED)
#SET(CUDA_VERBOSE_BUILD ON)

# ITK
#FIND_PACKAGE(ITK REQUIRED)
#INCLUDE(${ITK_USE_FILE})
#SET( ITK_LIBS ITKCommon ITKIO )

#INCLUDE_DIRECTORIES("../libtiff" ./cuda_common)
#LINK_DIRECTORIES("../libtiff")

#LINK_LIBRARIES( ${ITK_LIBS} )

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
    recon_api.h)

  SET_TARGET_PROPERTIES(Application_BackscatterReconstruction_Algorithm
    PROPERTIES
      COMPILE_DEFINITIONS USE_CUDA=1
  )
ELSE(CUDA_FOUND)
  ADD_LIBRARY(Application_BackscatterReconstruction_Algorithm # no main.cpp!
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
    recon_api.h)
ENDIF(CUDA_FOUND)

TARGET_LINK_LIBRARIES(Application_BackscatterReconstruction_Algorithm
  ${SCI_ITK_LIBRARIES}
)

