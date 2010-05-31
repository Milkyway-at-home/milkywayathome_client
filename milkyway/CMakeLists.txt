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

# The version number.
set(milkyway_client_VERSION_MAJOR 0)
set(milkyway_client_VERSION_MINOR 1)
set(milkyway_client_VERSION       0.1)

#if we aren't build from the root of the repo, these need to be set separately
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake_modules)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../bin)
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)

include_directories(${PROJECT_INCLUDE_DIR})

find_package(BOINC REQUIRED)
find_package(CUDA)

# if(APPLE)
#   find_library(OPENCL_LIBRARY OpenCL DOC "OpenCL lib for OS X")
#   find_path(OPENCL_INCLUDE_DIR OpenCL/cl.h DOC "Include for OpenCL on OS X")
# endif(APPLE)

include_directories(${BOINC_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) #FIXME: Move headers to include

include(CheckC99)
check_c99_support()

#CHECKME: Is this portable?
set(C99_FLAG -std=c99)

if(HAVE_C99_SUPPORT)
    set(CMAKE_C_FLAGS "${C99_FLAG} ${CMAKE_C_FLAGS}")
endif()


option(DOUBLE_PRECISION
       "Use double precision" ON)

option(MILKYWAY_GPU
       "Milkyway GPU" OFF)

option(COMPUTE_ON_GPU
       "Compute on GPU" OFF)

option(COMPUTE_ON_CPU
       "Compute on CPU" ON)

option(USE_OCL
       "Use OpenCL" OFF)

if(USE_OCL)
  include_directories ("${PROJECT_SOURCE_DIR}/milkyway_ocl")
endif(USE_OCL)

# Something about the preprocessor directives didn't seem right. It
# looks like an attempt at logic that belongs more in the build
# system. I think the intention was to turn on double precision by
# default for the CPU only, and have it be an option if you use the
# GPU so that's what I do. i.e. Double precision required for only
# CPU, optional otherwise.
if(COMPUTE_ON_CPU AND (NOT COMPUTE_ON_GPU))
  set(DOUBLE_PRECISION ON)
endif(COMPUTE_ON_CPU AND (NOT COMPUTE_ON_GPU))
#Fall through to using whatever the user set.


if(COMPUTE_ON_GPU)
    include_directories ("${PROJECT_SOURCE_DIR}/milkyway_gpu")
endif(COMPUTE_ON_GPU)

option(BOINC_APP_GRAPHICS
       "BOINC app graphics" OFF)



#FIXME: Somehow test for and require the Apple gcc.  The macports gccs
#seem to not work when you try to link against frameworks,
#CoreFoundation etc.
if(APPLE)
   find_library(CORE_FOUNDATION_LIBRARY CoreFoundation)
   set(OS_SPECIFIC_LIBS ${CORE_FOUNDATION_LIBRARY})
endif(APPLE)

set(source_file_list evaluation_state.c
                     simple_evaluator.c
                     ../util/io_util.c
                     stCoords.c
                     atSurveyGeometry.c
                     milkyway_complex.c
                     numericalIntegration.c
                     stVector.c
                     stMath.c
                     evaluation_optimized.c
                     parameters.c
                     star_points.c
                     astronomy_worker.c
                     milkyway_client_app.c)

set(extra_warnings "-Wshadow -Wpointer-arith -Wcast-align -common-Wwrite-strings -compare -fno-common")

set(CMAKE_C_FLAGS "-Wall -Wextra -pedantic ${extra_warnings} ${CMAKE_C_FLAGS}")

message("---------------------------------------------")
message("   Building:          ${PROJECT_NAME}")
message("   System:            ${CMAKE_SYSTEM}")
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
message("   ----")
message("   Build tool:        ${CMAKE_BUILD_TOOL}")
message("   C Compiler:        ${CMAKE_C_COMPILER}")
message("   CFLAGS:            ${CMAKE_C_FLAGS}")
message("   ----")
message("   Double precision:  ${DOUBLE_PRECISION}")
message("   OpenCL enabled:    ${USE_OCL}")
message("   Milkyway GPU:      ${MILKYWAY_GPU}")
message("   Compute on GPU:    ${COMPUTE_ON_GPU}")
message("   Compute on CPU:    ${COMPUTE_ON_CPU}")
message("   BOINC graphics:    ${BOINC_APP_GRAPHICS}")
message("---------------------------------------------")


configure_file (
  "${PROJECT_INCLUDE_DIR}/config.h.in"
  "${PROJECT_INCLUDE_DIR}/config.h"
  )


# add the executable
add_executable(milkyway_client ${source_file_list})
target_link_libraries(milkyway_client ${BOINC_LIBRARIES} ${OS_SPECIFIC_LIBS})

