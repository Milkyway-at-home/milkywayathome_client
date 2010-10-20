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

include(StaticBOINC)
include(CLKernelInline)
include(SpecialLibs)
include(BOINCBinNames)
include(PrintInfo)

set(separation_VERSION_MAJOR 0)
set(separation_VERSION_MINOR 52)
set(separation_VERSION "${separation_VERSION_MAJOR}.${separation_VERSION_MINOR}")


if(BOINC_RELEASE_NAME)
  get_boinc_bin_name("milkyway_separation" "${separation_VERSION}" "")
  set(separation_bin_name "${BOINC_BIN_NAME}")
else()
  set(separation_bin_name "milkyway_separation")
endif()


set(SEPARATION_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
include_directories(${SEPARATION_INCLUDE_DIR})


set(CMAKE_C_FLAGS "-std=gnu99 -warn-uninit-values ${CMAKE_C_FLAGS}")


option(SEPARATION_INLINE_KERNEL
       "Inline the CL kernel" OFF)
mark_as_advanced(SEPARATION_INLINE_KERNEL)

option(SMALL_TESTS
       "Run small tests" ON)

option(LARGE_TESTS
       "Run large tests" OFF)

option(AUX_LARGE_TESTS
       "Run large aux tests" OFF)

option(AUX_SMALL_TESTS
       "Run small aux tests" OFF)

maybe_static()

find_package(Libm REQUIRED)

if(BOINC_APPLICATION)
  find_package(BOINC)
  include_directories(${BOINC_INCLUDE_DIR})
endif()

if(SEPARATION_OPENCL)
  find_package(OpenCL)
  if(NOT OPENCL_FOUND)
    message(WARNING "Separation OpenCL was requested, but OpenCL not found. Disabling OpenCL.")
    set(SEPARATION_OPENCL OFF)
  else()
    include_directories(${OPENCL_DIRECTORIES})
  endif()
endif()

set_os_specific_libs(${SEPARATION_OPENCL})


if(NEED_POPT_DEPS)
    #These are only dependencies of popt needed for static linking
  find_package(Libintl REQUIRED)
  find_package(Libiconv REQUIRED)
endif()

find_package(Popt REQUIRED)
include_directories(${POPT_INCLUDE_DIR})

print_separator()
message("Separation build:")
print_libs()
print_separator()



configure_file (
  "${SEPARATION_INCLUDE_DIR}/separation_config.h.in"
  "${PROJECT_INCLUDE_DIR}/separation_config.h")

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
                         src/likelihood.c
                         src/calculated_constants.c
                         src/r_points.c)

set(separation_headers include/calculated_constants.h
                       include/integrals_common.h
                       include/separation_types.h
                       include/coordinates.h
                       include/separation.h
                       include/evaluation.h
                       include/io_util.h
                       include/star_points.h
                       include/evaluation_state.h
                       include/likelihood.h
                       include/gauss_legendre.h
                       include/parameters.h
                       include/separation_config.h.in
                       include/integrals.h
                       include/r_points.h
                       include/separation_constants.h)

set(separation_cl_headers include/setup_cl.h
                          include/run_cl.h
                          include/separation_cl_defs.h
                          include/separation_cl_buffers.h)



if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  list(APPEND lib_source_file_list ${cl_source_list})
  list(APPEND separation_headers ${separation_cl_headers})

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
else()
  list(APPEND lib_source_file_list src/integrals.c)
endif()

list(APPEND separation_link_libs ${OS_SPECIFIC_LIBS} milkyway ${LIBM_LIBRARY})

add_library(separation STATIC ${lib_source_file_list} ${separation_headers})
target_link_libraries(separation ${separation_link_libs})



set(exe_link_libs ${POPT_LIBRARY}
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

set(separation_main_src src/separation_main.c)

if(MSVC)
  list(APPEND separation_main_src src/main_win.cpp)
endif()

add_executable(${separation_bin_name} ${separation_main_src})
milkyway_link(${separation_bin_name} ${BOINC_APPLICATION} ${BOINC_RELEASE} "${exe_link_libs}")

add_subdirectory(tests EXCLUDE_FROM_ALL)


