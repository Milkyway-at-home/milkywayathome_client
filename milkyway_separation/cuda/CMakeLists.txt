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
cmake_minimum_required(VERSION 2.6)
project(milkyway_cuda)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_SHARED_LIBS
       "Set to OFF to build static libraries" ON)


# set(cuda_lib_src evaluation_gpu6_double.cu
#                  evaluation_gpu6_float.cu
#                  evaluation_gpu6_likelihood.cu
#                  evaluation_gpu6_double_const.cu
#                  evaluation_gpu6_double_global.cu
#                  evaluation_gpu6_double_shared.cu
#                  evaluation_gpu6.cu)

set(cuda_lib_src evaluation_gpu6.cu)

set(cuda_c_src gauss_legendre.c
               coords.c
               cpu_coords.c
               cpu_r_constants.c)


include_directories("${PROJECT_BINARY_DIR}")


#set(CUDA_VERBOSE_BUILD ON)

#set(CUDA_BUILD_EMULATION ON)

#set(CMAKE_C_FLAGS -m32)
#set(CMAKE_CXX_FLAGS -m32)

#set(CUDA_64_BIT_DEVICE_CODE OFF)    # Does not work OS X apparently for now

#set(CUDA_SDK_ROOT_DIR  "/Developer/GPU Computing/C/")
set(CUDA_SDK_ROOT_DIR  "/home/matt/NVIDIA_GPU_Computing_SDK/C/")

#set(CUDA_NVCC_FLAGS "-DDOUBLE_PRECISION=1 -DSHARED_MEMORY=1 -DGLOBAL_MEMORY=1")

#if (APPLE)
#set(CMAKE_OSX_ARCHITECTURES i386)
#set(CUDA_64_BIT_DEVICE_CODE OFF)
#endif (Apple)


#include_directories("/Developer/GPU Computing/C/common/inc")
#cuda_include_directories("/Developer/GPU Computing/C/common/inc")
cuda_include_directories("/home/matt/NVIDIA_GPU_Computing_SDK/C/common/inc")

cuda_add_library(milkyway_cuda ${cuda_lib_src} STATIC)

target_link_libraries(milkyway_cuda cutil cuda)

#cuda_add_executable(milkyway_cuda ${cuda_lib_src} "../milkyway/milkyway_client_app.c")


