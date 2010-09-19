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

include(StaticBOINC)
include(CLKernelInline)
include(SpecialLibs)

enable_testing()

set(separation_VERSION_MAJOR 0)
set(separation_VERSION_MINOR 5)
set(separation_VERSION "${separation_VERSION_MAJOR}.${separation_VERSION_MINOR}")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set_os_specific_libs()

include_directories(${PROJECT_INCLUDE_DIR})
set(CMAKE_C_FLAGS "-std=gnu99 -warn-uninit-values ${CMAKE_C_FLAGS}")


option(SEPARATION_INLINE_KERNEL
       "Inline the CL kernel" ON)


configure_file (
  "${PROJECT_INCLUDE_DIR}/separation_config.h.in"
  "${PROJECT_INCLUDE_DIR}/separation_config.h")

#if(UNIX AND NOT APPLE)  #Unix that isn't OS X
#    set(OS_SPECIFIC_LIBS pthread)
#endif()

set(cl_source_list src/separation_cl_buffers.c
                   src/setup_cl.c
                   src/run_cl.c
                   src/separation_cl_defs.c)

set(lib_source_file_list src/evaluation_state.c
                         src/io_util.c
                         src/gauss_legendre.c
                         src/evaluation.c
                         src/parameters.c
                         src/star_points.c
                         src/integrals.c
                         src/likelihood.c
                         src/calculated_constants.c
                         src/r_points.c)

if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  list(APPEND lib_source_file_list ${cl_source_list})

  if(SEPARATION_INLINE_KERNEL)
    message("Inlining kernel")
    set(c_kernel_dir "${PROJECT_SOURCE_DIR}/ckernels")
    inline_kernel("${PROJECT_SOURCE_DIR}/kernels/integrals.cl"
                  "integral_kernel"
                  "${c_kernel_dir}"
                  "${PROJECT_INCLUDE_DIR}" "${PROJECT_SOURCE_DIR}/src" "${MILKYWAY_INCLUDE_DIR}")
   include_directories(${c_kernel_dir})
   list(APPEND lib_source_file_list "${c_kernel_dir}/integral_kernel.c")
  endif()
endif()

set(separation_link_libs ${OS_SPECIFIC_LIBS} milkyway ${LIBM_LIBRARY})

add_library(separation STATIC ${lib_source_file_list})
target_link_libraries(separation ${separation_link_libs})



set(exe_link_libs ${LIBM_LIBRARY}
                  ${POPT_LIBRARY}
                  separation
                  ${separation_link_libs})

if(NEED_POPT_DEPS)
  # Popt dependencies for static linking For some reason, on OS X /
  # windows you need these dependencies but not Linux.
  list(APPEND exe_link_libs ${LIBICONV_LIBRARY}
                            ${LIBINTL_LIBRARY})
endif()

if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  list(APPEND exe_link_libs ${OPENCL_LIBRARIES})
endif()

set(separation_bin_name milkyway_separation)

add_executable(${separation_bin_name} src/milkyway_client_app.c)
milkyway_link(${separation_bin_name} ${BOINC_APPLICATION} ${BOINC_RELEASE} "${exe_link_libs}")

add_subdirectory(tests EXCLUDE_FROM_ALL)


