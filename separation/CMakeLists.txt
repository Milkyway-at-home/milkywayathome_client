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
include(CPUNameTest)

set(separation_VERSION_MAJOR 0)
set(separation_VERSION_MINOR 47)
set(separation_VERSION "${separation_VERSION_MAJOR}.${separation_VERSION_MINOR}")

get_info_from_processor_name()

if(BOINC_RELEASE_NAME)
  get_boinc_bin_name("milkyway_separation" "${separation_VERSION}" "")
  set(separation_bin_name "${BOINC_BIN_NAME}")
else()
  set(separation_bin_name "milkyway_separation")
endif()


set(SEPARATION_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
include_directories(${SEPARATION_INCLUDE_DIR})

if(SYSTEM_IS_PPC)
  set(CMAKE_C_FLAGS "-maltivec -mfused-madd -mhard-float -mabi=altivec ${CMAKE_C_FLAGS}")
endif()

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

option(SEPARATION_CRLIBM
      "Use crlibm double" OFF)
mark_as_advanced(SEPARATION_CRLIBM)

option(SEPARATION_FDLIBM
      "Use fdlibm double" OFF)
mark_as_advanced(SEPARATION_FDLIBM)

if(NOT DOUBLEPREC AND SEPARATION_CRLIBM)
  message(WARNING "Attempting to use single precision with crlibm. Disabling use of crlibm")
  set(SEPARTION_CRLIBM OFF)
endif()

if(NOT DOUBLEPREC AND SEPARATION_FDLIBM)
  message(WARNING "Attempting to use single precision with fdlibm. Disabling use of fdlibm")
  set(SEPARTION_FDLIBM OFF)
endif()


maybe_static()

if(NOT SEPARATION_FDLIBM)
  find_package(Libm REQUIRED)
  include_directories(${LIBM_INCLUDE_DIR})
  set(SEPARATION_LIBM_LIBRARY ${LIBM_LIBRARY})
else()
  # Replace standard libm with Fdlibm
  find_package(Fdlibm REQUIRED)
  include_directories(${FDLIBM_INCLUDE_DIR})
  set(SEPARATION_LIBM_LIBRARY ${FDLIBM_LIBRARY})
endif()

find_package(BOINC)

if(BOINC_APPLICATION)
  include_directories(${BOINC_INCLUDE_DIR})
endif()

find_package(OpenCL)

if(SEPARATION_OPENCL AND NOT OPENCL_FOUND)
  message(WARNING "Separation OpenCL was requested, but OpenCL not found. Disabling OpenCL.")
  set(SEPARATION_OPENCL OFF)
else()
  include_directories(${OPENCL_DIRECTORIES})
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
message("   Double precision:    ${DOUBLEPREC}")
message("   Separation crlibm:   ${SEPARATION_CRLIBM}")
print_libs()
print_separator()



configure_file (
  "${SEPARATION_INCLUDE_DIR}/separation_config.h.in"
  "${PROJECT_INCLUDE_DIR}/separation_config.h")

set(cl_source_list src/separation_cl_buffers.c
                   src/separation_binaries.c
                   src/setup_cl.c
                   src/run_cl.c)

set(lib_source_file_list src/evaluation_state.c
                         src/io_util.c
                         src/gauss_legendre.c
                         src/evaluation.c
                         src/parameters.c
                         src/star_points.c
                         src/likelihood.c
                         src/coordinates.c
                         src/integrals.c
                         src/calculated_constants.c
                         src/separation_utils.c
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
                       include/separation_utils.h
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
    set(SEPARATION_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
    set(SEPARATION_KERNEL_DIR ${PROJECT_SOURCE_DIR}/kernels)

    # Can figure out many of the header dependencies of kernel with C
    # preprocessor with -MM flag.  However, some of the files only purpose
    # is to include some other files in the kernel, so we have to exclude
    # those ones.
    set(all_cl_files
      ${MILKYWAY_INSTALL_INCLUDE_DIR}/milkyway_config.h
      ${MILKYWAY_INSTALL_INCLUDE_DIR}/separation_config.h

      ${MILKYWAY_INCLUDE_DIR}/milkyway_extra.h
      # This one not used in kernel
      #${MILKYWAY_INCLUDE_DIR}/milkyway_math_functions_standard.h

      ${MILKYWAY_INCLUDE_DIR}/milkyway_math_functions.h
      ${MILKYWAY_INCLUDE_DIR}/milkyway_math_functions_cl.h

      ${MILKYWAY_INCLUDE_DIR}/milkyway_math.h

      ${MILKYWAY_INCLUDE_DIR}/real.h
      ${MILKYWAY_INCLUDE_DIR}/milkyway_vector_types.h
      ${MILKYWAY_INCLUDE_DIR}/milkyway_vectors_func.h
      ${MILKYWAY_INCLUDE_DIR}/milkyway_vectors.h
      ${MILKYWAY_INCLUDE_DIR}/milkyway_cl.h

      ${SEPARATION_INCLUDE_DIR}/separation_types.h
      ${SEPARATION_INCLUDE_DIR}/integrals_common.h

      ${SEPARATION_KERNEL_DIR}/integrals.cl)

    set(INLINE_KERNEL_FILE "${LIBRARY_OUTPUT_PATH}/AllInOneFile.c")
    all_in_one_file("${all_cl_files}")
    inline_kernel("inlinedIntegralKernelSrc" "${ALL_IN_ONE_STRING}" ${INLINE_KERNEL_FILE})
    add_library(all_in_one_library STATIC ${INLINE_KERNEL_FILE})
    list(APPEND separation_link_libs all_in_one_library)
  endif()
endif()

list(APPEND separation_link_libs ${OS_SPECIFIC_LIBS} milkyway ${SEPARATION_LIBM_LIBRARY})

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

find_package(Crlibm)

if(CRLIBM_FOUND AND SEPARATION_CRLIBM)
  include_directories(${CRLIBM_INCLUDE_DIR})
  list(APPEND exe_link_libs ${CRLIBM_LIBRARY})
elseif(NOT CRLIBM_FOUND AND SEPARATION_CRLIBM)
  message(WARNING "Trying to build separation crlibm, but crlibm not found. Disabling crlibm")
  set(SEPARATION_CRLIBM OFF)
endif()

list(APPEND exe_link_libs ${SEPARATION_LIBM_LIBRARY})

if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  list(APPEND exe_link_libs ${OPENCL_LIBRARIES})
endif()

set(separation_main_src src/separation_main.c)

if(MSVC)
  list(APPEND separation_main_src src/main_win.cpp)
endif()

add_executable(${separation_bin_name} ${separation_main_src})
milkyway_link(${separation_bin_name} ${BOINC_APPLICATION} ${STATIC_EVERYTHING} "${exe_link_libs}")

add_subdirectory(tests EXCLUDE_FROM_ALL)




