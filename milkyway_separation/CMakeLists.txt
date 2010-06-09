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

# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

#No idea on the actual version, but I'll go with this for now.

cmake_minimum_required (VERSION 2.6)
project(milkyway_client)

include_directories("${PROJECT_SOURCE_DIR}")

option(BUILD_SHARED_LIBS
       "Set to OFF to build static libraries" ON)

option(DOUBLE_PRECISION
       "Use double precision" ON)

option(COMPUTE_ON_CPU
       "Compute on CPU" ON)

option(USE_CUDA
       "Use Nvidia CUDA" OFF)

option(USE_OCL
       "Use OpenCL" OFF)

#temporary stupidity to make sure you only pick one of the three.
#Build targets later.
if(COMPUTE_ON_CPU)
  list(APPEND appCount "CPU")
endif()

if(USE_CUDA)
  list(APPEND appCount "CUDA")
  include_directories("${CMAKE_SOURCE_DIR}/milkyway_separation_cuda")
  include_directories(${CUDA_INCLUDE_DIRS})
endif()

if(USE_OCL)
  list(APPEND appCount "OCL")
  include_directories("${CMAKE_SOURCE_DIR}/milkyway_separation_ocl")
endif()

list(LENGTH appCount NUM_CHOICES)

if(NOT (NUM_CHOICES EQUAL 1))
  message(FATAL_ERROR "Must choose one of CPU, OpenCL or CUDA.")
endif()

option(BOINC_APP_GRAPHICS
       "BOINC app graphics" OFF)


include(CheckAppleGCC)

#FIXME: Somehow test for and require the Apple gcc.  The macports gccs
#seem to not work when you try to link against frameworks,
#CoreFoundation etc. There's probably a better way to check and handle
#this.
if(APPLE)
  check_apple_gcc()
  if(NOT HAVE_APPLE_GCC)
    message(WARNING "Attempting to use a non-Apple GCC. Apple GCC is required.")
    message(FATAL_ERROR "You should specify a compiler with -DCMAKE_C_COMPILER=/path/to/compiler")
    # message(WARNING "Attempting to fall back to /usr/bin/gcc")
    # set(CMAKE_C_COMPILER /usr/bin/gcc)
    # message("LOLOLOL: ${CMAKE_C_COMPILER}")
    # check_apple_gcc()
    # if(NOT HAVE_APPLE_GCC)
    #   message(FATAL_ERROR "Fallback failed. Install an Apple compiler or manually specify one.")
    # endif()
  endif()

   find_library(CORE_FOUNDATION_LIBRARY CoreFoundation)
   set(OS_SPECIFIC_LIBS ${CORE_FOUNDATION_LIBRARY})
endif(APPLE)

if(UNIX AND NOT APPLE)  #Unix that isn't OS X
    set(OS_SPECIFIC_LIBS pthread)
endif(UNIX AND NOT APPLE)

if(WIN32)
  set(OS_SPECIFIC_LIBS msvcrt)
endif(WIN32)

set(lib_source_file_list evaluation_state.c
                         simple_evaluator.c
                         ../util/io_util.c
                         atSurveyGeometry.c
                         stCoords.c
                         stVector.c
                         numericalIntegration.c
                         evaluation_optimized.c
                         parameters.c
                         star_points.c
                         astronomy_worker.c)

set(gcc_extra_warnings "-Wshadow -Wpointer-arith -Wcast-align -common-Wwrite-strings -compare -fno-common")

#TODO: Other compilers
if(CMAKE_COMPILER_IS_GNUCC) # I think this could maybe be IS_UNIX
  set(CMAKE_C_FLAGS "-Wall -Wextra ${CMAKE_C_FLAGS}")
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
message("   BOINC Libraries    ${BOINC_LIBRARIES}")
message("   CUDA  Libraries    ${CUDA_LIBRARIES}")
message("   POPT  Library      ${POPT_LIBRARY}")
message("   ----")
message("   Build tool:        ${CMAKE_BUILD_TOOL}")
message("   C Compiler:        ${CMAKE_C_COMPILER}")
#TODO:Report CFLAGS used based on build type
#message("   CFLAGS:            ${CMAKE_C_FLAGS}")
message("   ----")
message("   Double precision:  ${DOUBLE_PRECISION}")
message("   OpenCL enabled:    ${USE_OCL}")
message("   CUDA enabled:      ${USE_CUDA}")
message("   Compute on CPU:    ${COMPUTE_ON_CPU}")
message("   BOINC graphics:    ${BOINC_APP_GRAPHICS}")
message("---------------------------------------------")


configure_file (
  "${PROJECT_INCLUDE_DIR}/config.h.in"
  "${PROJECT_INCLUDE_DIR}/config.h"
  )

# Windows is a pain in the ass.
if(WIN32)
  add_library(milkyway STATIC ${lib_source_file_list})
else(WIN32)
  add_library(milkyway SHARED ${lib_source_file_list})
endif(WIN32)

target_link_libraries(milkyway ${BOINC_LIBRARIES} ${OS_SPECIFIC_LIBS} m)

add_executable(milkyway_client milkyway_client_app.c)
target_link_libraries(milkyway_client ${POPT_LIBRARY} ${BOINC_LIBRARIES} milkyway)

install(TARGETS milkyway DESTINATION lib
       PERMISSIONS           #Fix occasional issuse with people who change their umask
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_INCLUDE_DIR}/milkyway.h DESTINATION include/milkyway)
install(FILES ${PROJECT_INCLUDE_DIR}/config.h DESTINATION include/milkyway)
install(TARGETS milkyway_client DESTINATION bin
       PERMISSIONS
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)

