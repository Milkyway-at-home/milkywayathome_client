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

#add_executable(kernel_inliner ${PROJECT_TOOL_DIR}/InlineKernel.c)


#file(READ zero_integral_kernel.cl zero_integral_kernel_str)
#file(READ zero_likelihood_kernel.cl zero_likelihood_kernel_str)

#file(READ nvidia_likelihood_kernel.cl nvidia_likelihood_kernel_str)

#message("HELLO: ${nvidia_likelihood_kernel_str}")

#"(/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/)|(//.*)"

#perl
#/\*(?:.|[\r\n])*?\*/

#string(REGEX REPLACE "(/\\*([^*]|[\r\n]|(\\*+([^*/]|[\r\n])))*\\*+/)|(//[^\r\n]*)" "" nvidia_likelihood_kernel_str_stripped ${nvidia_likelihood_kernel_str})

#file(WRITE "/tmp/lol.cl" ${nvidia_likelihood_kernel_str_stripped})

#string(REGEX REPLACE "(/\\*([^*]|[\r\n]|(\\*+([^*/]|[\r\n])))*\\*+/)|(//.*)" "" nvidia_likelihood_kernel_str_stripped ${nvidia_likelihood_kernel_str})

#string(REGEX REPLACE "(/\\*([^*]|[\r\n]|(\\*+([^*/]|[\r\n])))*\\*+/)|(//.*)" "" nvidia_likelihood_kernel_str_stripped ${nvidia_likelihood_kernel_str})

#message("HELLO: ${nvidia_likelihood_kernel_str_stripped}")
#message("-------------")

add_library(milkyway_separation_ocl STATIC ${ocl_lib_src})

target_link_libraries(milkyway_separation_ocl ${OPENCL_LIBRARIES})

# function(inline_kernel kernel_name)
#   add_custom_command(TARGET milkyway_separation_ocl
#                      PRE_BUILD
#                      ${kernel_inliner} "${kernel_name}.cl"
#                                        "${kernel_name}"
#                                        "${kernel_name}")
# endfunction()

# foreach(i ${kernel_list})
#   inline_kernel(${i})
# endforeach()


