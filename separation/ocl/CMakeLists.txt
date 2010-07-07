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

include(CLKernelInline)
#include_directories(${OPENCL_INCLUDE_DIRS})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories("../cpu/")
include_directories("../cuda/")


set(kernel_list ati_integral_kernel.cl
                ati_integral_kernel2.cl
                ati_likelihood_kernel.cl
                exp.cl
                nvidia_integral_kernel.cl
                nvidia_integral_kernel2.cl
                nvidia_likelihood_kernel.cl
                sqrt.cl
                zero_integral_kernel.cl
                zero_likelihood_kernel.cl)

set(ocl_lib_src evaluation_ocl.cpp
                evaluation_ocl_kernel.cpp
                evaluation_ocl_mem.cpp
                evaluation_ocl_setup.cpp
                ../cuda/gauss_legendre.c
                ../cuda/coords.c
                ../cuda/cpu_coords.c)

set(c_kernel_dir "${CMAKE_CURRENT_SOURCE_DIR}/ckernels")

#TODO: Bother only if the file has changed, and make sure it happens
#on every 'make'
foreach(kernel ${kernel_list})
  get_filename_component(basename "${kernel}" NAME_WE)
  inline_kernel("${basename}"
                "${CMAKE_CURRENT_SOURCE_DIR}/kernels/${kernel}"
                "${c_kernel_dir}")
  list(APPEND ocl_lib_src "${c_kernel_dir}/${basename}.c")
endforeach()

include_directories(${c_kernel_dir})

add_library(milkyway_separation_ocl STATIC ${ocl_lib_src})

target_link_libraries(milkyway_separation_ocl ${OPENCL_LIBRARIES})



