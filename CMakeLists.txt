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

cmake_minimum_required (VERSION 2.8)
project(milkywayathome)

enable_testing()


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

option(BOINC_APPLICATION
       "Run as BOINC application or independently" ON)

option(BOINC_RELEASE
       "Build with all options for making an actual release" ON)



#set(CMAKE_MODULE_PATH "${milkywayathome_SOURCE_DIR}/cmake_modules;${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "${milkywayathome_SOURCE_DIR}/cmake_modules")
set(PROJECT_TOOL_DIR "${milkywayathome_SOURCE_DIR}/tools")

#choose release build type by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MW_ENABLE_DEBUG ON)
  set(USE_BOINC_DIAGNOSTICS ON)
endif()

include_directories(${PROJECT_INCLUDE_DIR})

if(MINGW)
  if(NOT TAKEOFFGW_ROOT)
    set(TAKEOFFGW_ROOT "C:/takeoffgw/")
  endif(NOT TAKEOFFGW_ROOT)

  set(CMAKE_FIND_ROOT_PATH "${TAKEOFFGW_ROOT}/i686-pc-mingw32/sys-root/mingw")
endif(MINGW)

if(BOINC_RELEASE)
  set(POPT_USE_STATIC 1)
  set(JSON_C_USE_STATIC 1)
  set(LIBINTL_USE_STATIC 1)
  set(LIBICONV_USE_STATIC 1)
  set(BOINC_USE_STATIC 1)
endif()

if(UNIX AND NOT APPLE)
  set(NEED_POPT_DEPS 0)
else()
  set(NEED_POPT_DEPS 1)
endif()

#FIXME: relies on manually setting
if(NEED_POPT_DEPS)
    #These are only dependencies of popt needed for static linking
  find_package(Libintl REQUIRED)
  find_package(Libiconv REQUIRED)
endif()

list(APPEND CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_SOURCE_DIR}/nbody/lib")
list(APPEND CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/nbody/include")
find_package(Popt REQUIRED)
find_package(BOINC)
find_package(CUDA)
find_package(OpenCL)

#FIXME: Somehow test for and require the Apple gcc. The macports gccs
#seem to not work when you try to link against frameworks,
#CoreFoundation etc. There's probably a better way to check and handle
#this.
if(APPLE)
  include(CheckAppleGCC)
  check_apple_gcc()
  if(NOT HAVE_APPLE_GCC)
    message(WARNING "Attempting to use a non-Apple GCC. An Apple compiler is required on OS X.")
    message(FATAL_ERROR
             "You should specify a compiler by setting the CC environment variable or "
             "with -DCMAKE_C_COMPILER=/path/to/compiler. This should be /usr/bin/gcc usually.")
  endif()

  find_library(COREFOUNDATION_LIBRARY CoreFoundation )

  list(APPEND OS_SPECIFIC_LIBS ${COREFOUNDATION_LIBRARY})

  if(BOINC_RELEASE)
    # We don't want to bother compiling twice until we need the universal binary for releases
    #set(CMAKE_OSX_ARCHITECTURES i386;x86_64)
    set(CMAKE_OSX_ARCHITECTURES x86_64)
  endif()
endif(APPLE)

#Some kind of workaround for a popular cmake bug where make test
#doesn't behave as anyone would ever expect.
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

# If we found all the dependencies for something, by default try to
# build it


if(BOINC_FOUND AND POPT_FOUND)
  set(SEPARATION_IS_BUILDABLE ON)
  include_directories(${BOINC_INCLUDE_DIR})
  include_directories(${POPT_INCLUDE_DIR})

  if(CUDA_FOUND)
    set(CUDA_IS_BUILDABLE ON)
  endif()

  if(OPENCL_FOUND)
    set(OPENCL_IS_BUILDABLE ON)
    include_directories(${OPENCL_INCLUDE_DIRS})
  endif()
endif()

if(POPT_FOUND)
  set(NBODY_IS_BUILDABLE ON)
endif()


#You can turn things off if you really want to
option(CPU_SEPARATION
       "Build the cpu separation" ${SEPARATION_IS_BUILDABLE})

option(CUDA_SEPARATION
        "Build the CUDA separation" OFF)  # ${CUDA_IS_BUILDABLE})

option(OPENCL_SEPARATION
        "Build the OpenCL separation" ${OPENCL_IS_BUILDABLE})

option(NBODY
    "Build the n-body code" ${NBODY_IS_BUILDABLE})

option(DOUBLE_PRECISION
         "Use double precision" ON)

option(BOINC_APP_GRAPHICS
       "BOINC app graphics" OFF)


set(gcc_extra_warnings "-Wshadow -Wpointer-arith -Wcast-align -common-Wwrite-strings -compare -fno-common -Wstrict-aliasing -Wstrict-aliasing=3 -Wswitch-enum -Wswitch-default -Wmissing-declarations -Wfloat-equal -Wcomment")

#TODO: Other compilers
if(CMAKE_COMPILER_IS_GNUCC) # I think this could maybe be IS_UNIX
  set(CMAKE_C_FLAGS "-Wall -Wextra ${CMAKE_C_FLAGS} -funroll-loops -maccumulate-outgoing-args")

  set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra ${gcc_extra_warnings} -O0 -g ${CMAKE_C_FLAGS_DEBUG}")
endif(CMAKE_COMPILER_IS_GNUCC)



message("---------------------------------------------")
message("   Building:          ${PROJECT_NAME}")
message("   System:            ${CMAKE_SYSTEM}")
message("   Build type:        ${CMAKE_BUILD_TYPE}")
message("   Arch:              ${CMAKE_SYSTEM_PROCESSOR}")
message("   Install path:      ${CMAKE_INSTALL_PREFIX}")
message("   ----")
message("   CMAKE version:     ${CMAKE_VERSION}")
message("   CMAKE binary:      ${CMAKE_COMMAND}")
message("   CTEST binary:      ${CMAKE_CTEST_COMMAND}")
message("   CMAKE generator:   ${CMAKE_GENERATOR}")
message("   ----")
message("   Project src dir:   ${CMAKE_SOURCE_DIR}")
message("   Project bin dir:   ${CMAKE_BINARY_DIR}")
message("   ----")
message("   BOINC Library:     ${BOINC_LIBRARY}")
message("   CUDA  Libraries    ${CUDA_LIBRARIES}")
message("   OpenCL Libraries   ${OPENCL_LIBRARIES}")
message("   POPT Library       ${POPT_LIBRARY}")
message("   ----")
message("   Build tool:        ${CMAKE_BUILD_TOOL}")
message("   C Compiler:        ${CMAKE_C_COMPILER}")
#TODO:Report CFLAGS used based on build type
#message("   CFLAGS:           ${CMAKE_C_FLAGS}")
message("   ----")
message("   Double precision:  ${DOUBLE_PRECISION}")
message("   OpenCL separation: ${OPENCL_SEPARATION}")
message("   CUDA separation:   ${CUDA_SEPARATION}")
message("   CPU separaton:     ${CPU_SEPARATION}")
message("   N-Body:            ${NBODY}")
message("   BOINC graphics:    ${BOINC_APP_GRAPHICS}")
message("---------------------------------------------")


if(NOT CPU_SEPARATION AND NOT CUDA_SEPARATION AND NOT OPENCL_SEPARATION AND NOT NBODY)
  message(FATAL_ERROR "Nothing will be built!")
endif()

if(SEPARATION_IS_BUILDABLE)
  add_subdirectory(separation)
endif()

if(NBODY_IS_BUILDABLE)
  add_subdirectory(nbody)
endif()


