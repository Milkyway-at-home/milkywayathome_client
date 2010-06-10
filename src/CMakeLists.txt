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

project(milkyway_client)


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
  include_directories("${CMAKE_SOURCE_DIR}/milkyway_separation_cpu")
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


include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_SOURCE_DIR}/milkyway_separation/cpu")


set(USE_CPU_DEFINES "-DCOMPUTE_ON_CPU=1 -DUSE_CUDA=0 -D USE_OCL=0")


# CPU Variant
add_executable(milkyway_client milkyway_client_app.c)
set_target_properties(milkyway_client
                        PROPERTIES
                          COMPILE_FLAGS "-DCOMPUTE_ON_CPU=1 -DUSE_CUDA=0 -D USE_OCL=0")


target_link_libraries(milkyway_client ${POPT_LIBRARY} ${BOINC_LIBRARIES} milkyway_separation_cpu)

install(TARGETS milkyway_client DESTINATION bin
       PERMISSIONS
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)


# OpenCL variant
add_executable(milkyway_client_ocl milkyway_client_app.c)
target_link_libraries(milkyway_client_ocl ${POPT_LIBRARY}
                                          ${BOINC_LIBRARIES}
                                          milkyway_separation_ocl
                                          milkyway_separation_cpu)
set_target_properties(milkyway_client_ocl
                        PROPERTIES
                          COMPILE_FLAGS "-DCOMPUTE_ON_CPU=0 -DUSE_CUDA=0 -D USE_OCL=1")



install(TARGETS milkyway_client_ocl DESTINATION bin
       PERMISSIONS
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)

