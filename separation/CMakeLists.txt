# Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
# Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
# Magdon-Ismail and Rensselaer Polytechnic Institute.
#
# This file is part of Milkway@Home.
#
# Milkyway@Home is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Milkyway@Home is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.6)
project(milkyway_separation)

enable_testing()

include(StaticBOINC)
include(CLKernelInline)
include(SpecialLibs)
include(BOINCBinNames)
include(PrintInfo)
include(CPUNameTest)
include(BuildSSELevels)
include(CMakeDependentOption)
include(BOINCPaths)

set(separation_VERSION_MAJOR 0
               CACHE INTERNAL "Separation version number")
set(separation_VERSION_MINOR 88
               CACHE INTERNAL "Separation version number")
set(separation_VERSION "${separation_VERSION_MAJOR}.${separation_VERSION_MINOR}"
               CACHE INTERNAL "Separation version number")

include_directories(${MILKYWAY_INCLUDE_DIR})
include_directories(${MILKYWAY_INSTALL_INCLUDE_DIR})
include_directories(${DSFMT_INCLUDE_DIR})

set(SEPARATION_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "Sepaation includes")
set(SEPARATION_KERNEL_DIR "${PROJECT_SOURCE_DIR}/kernels" CACHE INTERNAL "Separation kernels")
include_directories(${SEPARATION_INCLUDE_DIR})

if(SYSTEM_IS_PPC)
  set(CMAKE_C_FLAGS "-maltivec -mfused-madd -mhard-float -mabi=altivec ${CMAKE_C_FLAGS}")
endif()

option(SEPARATION_STATIC "Build separation as fully static binary" ON)
mark_as_advanced(SEPARATION_STATIC)

option(SEPARATION_CAL "Use CAL++ voodoo" OFF)

#FIXME: Still need to make it re-inline if kernel changed
cmake_dependent_option(SEPARATION_INLINE_KERNEL "Inline the CL kernel" ON
                         "SEPARATION_OPENCL" OFF)
mark_as_advanced(SEPARATION_INLINE_KERNEL)

option(SEPARATION_GRAPHICS
       "Use Separation Graphics" OFF)

cmake_dependent_option(SEPARATION_CRLIBM "Use crlibm for separation maths" OFF
                       "DOUBLEPREC" OFF)
mark_as_advanced(SEPARATION_CRLIBM)

cmake_dependent_option(SEPARATION_FDLIBM "Use fdlibm for separation maths" OFF
                       "DOUBLEPREC" OFF)
mark_as_advanced(SEPARATION_FDLIBM)

if(NOT DOUBLEPREC AND SEPARATION_CRLIBM)
  message(WARNING "Attempting to use single precision with crlibm. Disabling use of crlibm")
  set(SEPARTION_CRLIBM OFF)
endif()

if(NOT DOUBLEPREC AND SEPARATION_FDLIBM)
  message(WARNING "Attempting to use single precision with fdlibm. Disabling use of fdlibm")
  set(SEPARTION_FDLIBM OFF)
endif()

if(NOT BOINC_APPLICATION AND SEPARATION_GRAPHICS)
  message(WARNING "Attempting to separation graphics without BOINC. Disabling graphics")
  set(SEPARATION_GRAPHICS OFF)
endif()

if(SEPARATION_CAL AND (NOT CAL++_FOUND OR NOT CAL_FOUND))
  message(FATAL_ERROR "Separation CAL was requested, but CAL or CAL++ not found.")
endif()

if(SEPARATION_OPENCL)
  message(FATAL_ERROR "Separation OpenCL was requested, but OpenCL not found.")
else()
  include_directories(${OPENCL_DIRECTORIES})
endif()

if(SEPARATION_OPENCL AND SEPARATION_CAL)
  message(FATAL_ERROR "Can only choose one of OpenCL or CAL at a time")
endif()


if(SEPARATION_CAL OR SEPARATION_OPENCL)
  set(SEPARATION_GPU 1)
endif()

if(SEPARATION_GPU)
  message(AUTHOR_WARNING "Using GPU, so unsetting static separation")
  set(SEPARATION_STATIC 0)
endif()

maybe_static(${SEPARATION_STATIC})


if(NOT FDLIBM_FOUND AND SEPARATION_FDLIBM)
  message(WARNING "Separation Fdlibm was requested, but Fdlibm found. Disabling Fdlibm.")
  set(SEPARATION_FDLIBM OFF)
endif()

if(NOT SEPARATION_FDLIBM)
  if(NOT WIN32)
    set(SEPARATION_LIBM_LIBRARY m)
  endif()
else()
  # Replace standard libm with Fdlibm
  include_directories(${FDLIBM_INCLUDE_DIR})
  set(SEPARATION_LIBM_LIBRARY ${FDLIBM_LIBRARY})
endif()



if(BOINC_APPLICATION)
  include_directories(${BOINC_INCLUDE_DIR})
endif()

include_directories(${POPT_INCLUDE_DIR})


print_separator()
message("Separation build:")
message("   Double precision:    ${DOUBLEPREC}")
message("   Separation crlibm:   ${SEPARATION_CRLIBM}")
message("   Separation OpenCL:   ${SEPARATION_OPENCL}")
message("   Separation CAL:      ${SEPARATION_CAL}")
print_libs()
print_separator()

configure_file (
  "${SEPARATION_INCLUDE_DIR}/separation_config.h.in"
  "${PROJECT_INCLUDE_DIR}/separation_config.h")

set(cl_source_list src/separation_cl_buffers.c
                   src/separation_binaries.c
                   src/setup_cl.c
                   src/run_cl.c)

set(cal_source_list #src/separation_cal.cpp
                    src/separation_cal_kernelgen.cpp
                    src/separation_cal_setup.c
                    src/separation_cal_run.c
                    src/show_cal_types.c)


set(graphics_src_list src/separation_graphics.cc)
set(graphics_hdr_list include/separation_graphics.h)

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
                         src/r_points.c
                         src/separation_lua.c)

set(separation_headers include/calculated_constants.h
                       include/separation_types.h
                       include/separation_kernel_types.h
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
                       include/separation_constants.h
                       include/separation_lua.h)

set(separation_cl_headers include/setup_cl.h
                          include/run_cl.h
                          include/separation_cl_defs.h
                          include/separation_cl_buffers.h)

set(separation_cal_headers include/separation_cal.h
                           include/separation_cal_kernelgen.h
                           include/separation_cal_types.h
                           include/separation_cal_setup.h
                           include/separation_cal_run.h
                           include/show_cal_types.h)

if(SEPARATION_GRAPHICS)
  list(APPEND lib_source_file_list ${graphics_src_list})
  list(APPEND separation_headers ${graphics_hdr_list})
endif()

if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  list(APPEND lib_source_file_list ${cl_source_list})
  list(APPEND separation_headers ${separation_cl_headers})

  if(SEPARATION_INLINE_KERNEL)
    # Can figure out many of the header dependencies of kernel with C
    # preprocessor with -MM flag.  However, some of the files only purpose
    # is to include some other files in the kernel, so we have to exclude
    # those ones.
    set(all_cl_files
          ${SEPARATION_INCLUDE_DIR}/separation_kernel_types.h
          ${SEPARATION_KERNEL_DIR}/integrals.cl)


    set(SEPARATION_INLINE_KERNEL_FILE "${LIBRARY_OUTPUT_PATH}/AllInOneFile.c")
    set(SEPARATION_INLINE_SCRIPT "${CMAKE_MODULE_PATH}/InlineSeparation.cmake")
    set(SEPARATION_OUT_STR_NAME "inlinedIntegralKernelSrc")
    add_custom_command(
      OUTPUT "${SEPARATION_INLINE_KERNEL_FILE}"
      COMMAND ${CMAKE_COMMAND} "-DOUTFILE=${SEPARATION_INLINE_KERNEL_FILE}"
                               "-DSOURCEFILES=\"${all_cl_files}\""
                               "-DSTRNAME=${SEPARATION_OUT_STR_NAME}"
                               "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}"
                               -P ${SEPARATION_INLINE_SCRIPT})
    add_library(all_in_one_library STATIC "${SEPARATION_INLINE_KERNEL_FILE}")
    list(APPEND separation_link_libs all_in_one_library)
  endif()
elseif(SEPARATION_CAL)
  include_directories(${Boost_INCLUDE_DIR})

  list(APPEND lib_source_file_list ${cal_source_list})
  list(APPEND separation_headers ${separation_cal_headers})
endif()

list(APPEND separation_link_libs milkyway milkyway_lua)

if(BOINC_APPLICATION)
  list(APPEND separation_link_libs ${BOINC_LIBRARIES})
endif()

list(APPEND separation_link_libs ${OS_SPECIFIC_LIBS}
                                 ${LUA_LIBRARIES}
                                 ${SEPARATION_LIBM_LIBRARY})


set(exe_link_libs ${POPT_LIBRARY}
                  ${separation_link_libs})

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
  list(APPEND exe_link_libs milkyway_cl ${OPENCL_LIBRARIES})
elseif(SEPARATION_CAL)
  include_directories(${ATICAL_INCLUDE_DIR})
  include_directories(${ATICAL_INCLUDE_DIR}/CAL)
  include_directories(${CAL++_INCLUDE_DIR})
  include_directories(${Boost_INCLUDE_DIR})
  list(APPEND exe_link_libs ${ATICALCL_LIBRARY} ${ATICALRT_LIBRARY})
endif()

set(separation_main_src src/separation_main.c)

if(MSVC)
  list(APPEND separation_main_src src/main_win.c)
endif()

if(SEPARATION_OPENCL)
  if(NVIDIA_OPENCL)
    set(CL_PLAN_CLASS "cuda_opencl")
  elseif(AMD_OPENCL)
    set(CL_PLAN_CLASS "ati_opencl")
  else()
    set(CL_PLAN_CLASS "opencl")
  endif()
elseif(SEPARATION_CAL)
  set(CL_PLAN_CLASS "ati14")
else()
  set(CL_PLAN_CLASS "")
endif()

add_library(separation STATIC ${lib_source_file_list} ${separation_headers})
maybe_disable_ssen(separation)



if(NOT ALWAYS_HAVE_SSE2)
  add_library(separation_core STATIC src/probabilities.c include/probabilities.h)
  if(SYSTEM_IS_X86)
    disable_sse2(separation_core)
    disable_sse3(separation_core)
  endif()
  set(separation_core_libs separation_core)
endif()


# Build the critical area separately for SSE2, SSE3, no SSE2
if(SYSTEM_IS_X86)
  add_library(separation_core_sse2 STATIC src/probabilities.c include/probabilities.h)
  list(APPEND separation_core_libs separation_core_sse2)
  enable_sse2(separation_core_sse2)
  disable_sse3(separation_core_sse2)

  add_library(separation_core_sse3 STATIC src/probabilities.c include/probabilities.h)
  enable_sse2(separation_core_sse3)
  enable_sse3(separation_core_sse3)
  list(APPEND separation_core_libs separation_core_sse3)
endif()


if(NOT MILKYWAY_IPHONE_APP)
  add_executable(milkyway_separation ${separation_main_src})
  maybe_disable_ssen(milkyway_separation)
else()
  add_library(milkyway_separation STATIC ${separation_main_src})
endif()

if(BOINC_RELEASE_NAMES)
  get_boinc_bin_name("milkyway_separation" "${separation_VERSION}" "${CL_PLAN_CLASS}")
  set_target_properties(milkyway_separation
                          PROPERTIES
                            OUTPUT_NAME "${BOINC_BIN_NAME}")
endif()

milkyway_link(milkyway_separation ${BOINC_APPLICATION}
                                  ${SEPARATION_STATIC}
                                  "separation;${separation_core_libs};${exe_link_libs}")


add_subdirectory(tests EXCLUDE_FROM_ALL)

install(TARGETS milkyway_separation
          RUNTIME DESTINATION ${MILKYWAY_PROJECT_DIRECTORY})


