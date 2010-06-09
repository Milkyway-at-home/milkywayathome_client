# Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
# Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
# Magdon-Ismail and Rensselaer Polytechnic Institute.

# This file is part of Milkway@Home.

# Milkyway@Home is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Milkyway@Home is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required (VERSION 2.6)
project(milkyway@home)

enable_testing()

set(milkyway_client_VERSION_MAJOR 0)
set(milkyway_client_VERSION_MINOR 1)
set(milkyway_client_VERSION       0.1)


#choose release build type by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MW_ENABLE_DEBUG ON)
endif()


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

include_directories(${PROJECT_INCLUDE_DIR})

find_package(CUDA)

if(MINGW)
  if(NOT TAKEOFFGW_ROOT)
    set(TAKEOFFGW_ROOT "C:/takeoffgw/")
  endif(NOT TAKEOFFGW_ROOT)

  set(CMAKE_FIND_ROOT_PATH "${TAKEOFFGW_ROOT}/i686-pc-mingw32/sys-root/mingw")
endif(MINGW)

find_package(Popt REQUIRED)
find_package(BOINC REQUIRED)


# if(APPLE)
#   find_library(OPENCL_LIBRARY OpenCL DOC "OpenCL lib for OS X")
#   find_path(OPENCL_INCLUDE_DIR OpenCL/cl.h DOC "Include for OpenCL on OS X")
# endif(APPLE)

include_directories(${BOINC_INCLUDE_DIR})
include_directories(${POPT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) #FIXME: Move headers to include
include_directories(${PROJECT_INCLUDE_DIR})


#Some kind of workaround for a popular cmake bug where make test
#doesn't behave as anyone would ever expect.
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

option(SMALL_TESTS
       "Run small tests" ON)

option(LARGE_TESTS
       "Run large tests" OFF)

option(AUX_LARGE_TESTS
       "Run large aux tests" OFF)

option(AUX_SMALL_TESTS
       "Run small aux tests" OFF)

option(MW_DEBUG_EXCESSIVE
       "Excessive debugging output" OFF)


add_subdirectory(milkyway_separation)

#if(CUDA_FOUND)
#  add_subdirectory(milkyway_separation_cuda)
#endif()

if(OCL_FOUND)
  add_subdirectory(milkyway_separation_ocl)
endif()


add_subdirectory(tests EXCLUDE_FROM_ALL)


