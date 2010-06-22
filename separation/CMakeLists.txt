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

cmake_minimum_required(VERSION 2.8)
project(milkyway_separation)

enable_testing()

set(milkyway_separation_VERSION_MAJOR 0)
set(milkyway_separation_VERSION_MINOR 4)
set(milkyway_separation_${milkyway_client_MAJOR_VERSION}.${milkyway_client_MINOR_VERSION})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include_directories(${PROJECT_INCLUDE_DIR})

configure_file (
  "${PROJECT_INCLUDE_DIR}/config.h.in"
  "${PROJECT_INCLUDE_DIR}/config.h")


if(UNIX AND NOT APPLE)  #Unix that isn't OS X
    set(OS_SPECIFIC_LIBS pthread)
endif()

if(WIN32)
  set(OS_SPECIFIC_LIBS msvcrt)
endif()


if(CPU_SEPARATION)
  add_subdirectory(cpu)
endif()

if(CUDA_SEPARATION)
  add_subdirectory(cuda)
endif()

if(OPENCL_SEPARATION)
  add_subdirectory(ocl)
endif()

add_subdirectory(src)
add_subdirectory(tests EXCLUDE_FROM_ALL)


