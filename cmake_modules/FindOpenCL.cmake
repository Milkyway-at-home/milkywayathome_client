# Taken from http://gitorious.org/findopencl/findopencl 6/15/2010
# Says it's public domain

# - Try to find OpenCL
# This module tries to find an OpenCL implementation on your system. It supports
# AMD / ATI, Apple and NVIDIA implementations, but shoudl work, too.
#
# Once done this will define
#  OPENCL_FOUND        - system has OpenCL
#  OPENCL_INCLUDE_DIRS  - the OpenCL include directory
#  OPENCL_LIBRARIES    - link these to use OpenCL
#

include(CPUNameTest)

find_package( PackageHandleStandardArgs )

set(OPENCL_VERSION_STRING "0.1.0")
set(OPENCL_VERSION_MAJOR 0)
set(OPENCL_VERSION_MINOR 1)
set(OPENCL_VERSION_PATCH 0)

if(APPLE)

  FIND_LIBRARY(OPENCL_LIBRARIES OpenCL DOC "OpenCL lib for OSX")
  FIND_PATH(OPENCL_INCLUDE_DIRS OpenCL/cl.h DOC "Include for OpenCL on OSX")
  FIND_PATH(_OPENCL_CPP_INCLUDE_DIRS OpenCL/cl.hpp DOC "Include for OpenCL CPP bindings on OSX")

else(APPLE)

  if(WIN32)
    if(NVIDIA_OPENCL)
      # We could use CUDA_LIB_PATH, but this breaks when compiling 32 on 64
      if(SYSTEM_IS_64)
        set(OPENCL_LIB_DIR "$ENV{CUDA_PATH}/lib/x64/")
      else()
        set(OPENCL_LIB_DIR "$ENV{CUDA_PATH}/lib/Win32")
      endif()
      set(_OPENCL_INC_CAND "$ENV{CUDA_INC_PATH}")
    elseif(AMD_OPENCL)
      # The AMD SDK currently installs both x86 and x86_64 libraries
      if(SYSTEM_IS_64)
        set(OPENCL_LIB_DIR "$ENV{AMDAPPSDKROOT}/lib/x86_64")
      else()
        set(OPENCL_LIB_DIR "$ENV{AMDAPPSDKROOT}/lib/x86")
      endif()
      get_filename_component(_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include ABSOLUTE)
    endif()

    if(SYSTEM_IS_64)
      find_library(OPENCL_LIBRARIES OpenCL.lib ${OPENCL_LIB_DIR} "$ENV{CUDA_PATH}/lib/x64/" "$ENV{AMDAPPSDKROOT}/lib/x86_64")
    else()
      find_library(OPENCL_LIBRARIES OpenCL.lib ${OPENCL_LIB_DIR} "$ENV{CUDA_PATH}/lib/Win32" "$ENV{AMDAPPSDKROOT}/lib/x86")
    endif()

    # On Win32 search relative to the library
    # If we are forcing AMD or Nvidia OpenCL, those should be found before trying both
    find_path(OPENCL_INCLUDE_DIRS CL/cl.h PATHS "${_OPENCL_INC_CAND}" "$ENV{AMDAPPSDKROOT}/include" "$ENV{CUDA_INC_PATH}")
    find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS "${_OPENCL_INC_CAND}")
  else(WIN32)
    if(NVIDIA_OPENCL)
      find_library(OPENCL_LIBRARIES OpenCL
                  /usr/local/cuda)
    elseif(AMD_OPENCL)
      find_library(OPENCL_LIBRARIES OpenCL
                    $ENV{AMDAPPSDKROOT}/lib/x86_64
                    $ENV{AMDAPPSDKROOT}/lib/x86)
    else()
      find_library(OPENCL_LIBRARIES OpenCL
        $ENV{AMDAPPSDKROOT}/lib/x86_64
        $ENV{AMDAPPSDKROOT}/lib/x86
        /usr/local/cuda)
    endif()
    get_filename_component(OPENCL_LIB_DIR ${OPENCL_LIBRARIES} PATH)
    get_filename_component(_OPENCL_INC_CAND ${OPENCL_LIB_DIR}/../../include ABSOLUTE)

    # The AMD SDK currently does not place its headers
    # in /usr/include, therefore also search relative
    # to the library
    find_path(OPENCL_INCLUDE_DIRS CL/cl.h
                PATHS
                  ${_OPENCL_INC_CAND}
                  $ENV{AMDAPPSDKROOT}/include
                  /usr/local/cuda/include)
    find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp
                PATHS
                  ${_OPENCL_INC_CAND})
  endif (WIN32)

endif(APPLE)

find_package_handle_standard_args( OpenCL DEFAULT_MSG OPENCL_LIBRARIES OPENCL_INCLUDE_DIRS )

if( _OPENCL_CPP_INCLUDE_DIRS )
  set( OPENCL_HAS_CPP_BINDINGS TRUE )
  list( APPEND OPENCL_INCLUDE_DIRS ${_OPENCL_CPP_INCLUDE_DIRS} )
  # This is often the same, so clean up
  list( REMOVE_DUPLICATES OPENCL_INCLUDE_DIRS )
endif( _OPENCL_CPP_INCLUDE_DIRS )

mark_as_advanced(
  OPENCL_INCLUDE_DIRS
  )

if(OPENCL_INCLUDE_DIRS AND OPENCL_LIBRARIES)
  set(OPENCL_FOUND TRUE)
else()
  set(OPENCL_FOUND FALSE)
endif()


