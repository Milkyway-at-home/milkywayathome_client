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

include(MWLink)
include(SpecialLibs)
include(BOINCBinNames)
include(PrintInfo)
include(CPUNameTest)
include(CheckSSE2)
include(CMakeDependentOption)
include(BOINCPaths)

set(separation_VERSION_MAJOR 1
               CACHE INTERNAL "Separation version number")
set(separation_VERSION_MINOR 44
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

option(SEPARATION_STATIC "Build separation as fully static binary" OFF)
mark_as_advanced(SEPARATION_STATIC)



option(SEPARATION_GRAPHICS
       "Use Separation Graphics" OFF)

cmake_dependent_option(SEPARATION_CRLIBM "Use crlibm for separation maths" OFF
                       "DOUBLEPREC" OFF)
mark_as_advanced(SEPARATION_CRLIBM)

if(NOT DOUBLEPREC AND SEPARATION_CRLIBM)
  message(WARNING "Attempting to use single precision with crlibm. Disabling use of crlibm")
  set(SEPARTION_CRLIBM OFF)
endif()

if(NOT BOINC_APPLICATION AND SEPARATION_GRAPHICS)
  message(WARNING "Attempting to separation graphics without BOINC. Disabling graphics")
  set(SEPARATION_GRAPHICS OFF)
endif()

if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()


cmake_dependent_option(SEPARATION_STATIC "Build separation as fully static binary" OFF
                                         "NOT SEPARATION_OPENCL" OFF)
mark_as_advanced(SEPARATION_STATIC)
maybe_static(${SEPARATION_STATIC})

cmake_dependent_option(MSVC32_AVX_WORKAROUND
                       "Disable AVX build in case trying to build for Windows x86 with VS2008" ON
                       "SYSTEM_IS_X86_32;NOT SEPARATION_OPENCL;MSVC" OFF)

if(NOT WIN32)
  set(SEPARATION_LIBM_LIBRARY m)
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
print_libs()
print_separator()

configure_file (
  "${SEPARATION_INCLUDE_DIR}/separation_config.h.in"
  "${PROJECT_INCLUDE_DIR}/separation_config.h")

set(cl_source_list src/separation_cl_buffers.c
                   src/separation_binaries.c
                   src/cl_compile_flags.cpp
                   src/replace_amd_il.c
                   src/setup_cl.c
                   src/run_cl.c
                   src/il_kernel_1.c
                   src/il_kernel_2.c
                   src/il_kernel_3.c
                   src/il_kernel_4.c)

set(graphics_src_list src/separation_graphics.cc)
set(graphics_hdr_list include/separation_graphics.h)

set(lib_source_file_list src/evaluation_state.c
                         src/io_util.c
                         src/gauss_legendre.c
                         src/evaluation.c
                         src/probabilities.c
                         src/probabilities_dispatch.c
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
                       include/coordinates.h
                       include/separation.h
                       include/evaluation.h
                       include/probabilities_dispatch.h
                       include/probabilities.h
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
                          include/cl_compile_flags.h
                          include/replace_amd_il.h
                          include/il_kernels.h
                          include/run_cl.h
                          include/separation_cl_buffers.h)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/../boinc/boinc)

if(SEPARATION_GRAPHICS)
  list(APPEND lib_source_file_list ${graphics_src_list})
  list(APPEND separation_headers ${graphics_hdr_list})
endif()

if(SEPARATION_OPENCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
  list(APPEND lib_source_file_list ${cl_source_list})
  list(APPEND separation_headers ${separation_cl_headers})
  set(separation_inline_src "${LIBRARY_OUTPUT_PATH}/separation_kernel_inline.c")
  if(CrossCompile)
    add_custom_command(
      OUTPUT "${separation_inline_src}"
      COMMAND xxd -i "probabilities_kernel.cl" ${separation_inline_src}
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/kernels"
      DEPENDS "${PROJECT_SOURCE_DIR}/kernels/probabilities_kernel.cl"
      MAIN_DEPENDENCY "${PROJECT_SOURCE_DIR}/kernels/probabilities_kernel.cl"
      COMMENT "Inlining separation kernel source")
  else()
    add_custom_command(
      OUTPUT "${separation_inline_src}"
      COMMAND $<TARGET_FILE:xxd> -i "probabilities_kernel.cl" ${separation_inline_src}
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/kernels"
      DEPENDS "${PROJECT_SOURCE_DIR}/kernels/probabilities_kernel.cl"
      MAIN_DEPENDENCY "${PROJECT_SOURCE_DIR}/kernels/probabilities_kernel.cl"
      COMMENT "Inlining separation kernel source")
  endif()

  set(summarization_inline_src "${LIBRARY_OUTPUT_PATH}/summarization_kernel_inline.c")
  if(CrossCompile)
    add_custom_command(
      OUTPUT "${summarization_inline_src}"
      COMMAND xxd -i "summarization_kernel.cl" ${summarization_inline_src}
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/kernels"
      DEPENDS "${PROJECT_SOURCE_DIR}/kernels/summarization_kernel.cl"
      MAIN_DEPENDENCY "${PROJECT_SOURCE_DIR}/kernels/summarization_kernel.cl"
      COMMENT "Inlining summarization kernel source")
  else()
    add_custom_command(
      OUTPUT "${summarization_inline_src}"
      COMMAND $<TARGET_FILE:xxd> -i "summarization_kernel.cl" ${summarization_inline_src}
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/kernels"
      DEPENDS "${PROJECT_SOURCE_DIR}/kernels/summarization_kernel.cl"
      MAIN_DEPENDENCY "${PROJECT_SOURCE_DIR}/kernels/summarization_kernel.cl"
      COMMENT "Inlining summarization kernel source")
  endif()

  set(kernel_srcs ${separation_inline_src} ${summarization_inline_src})
  add_library(separation_kernels STATIC ${kernel_srcs})
  list(APPEND separation_link_libs separation_kernels)
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
  include_directories(${LIBELF_INCLUDE_DIR})
  list(APPEND exe_link_libs milkyway_cl ${OPENCL_LIBRARIES} ${LIBELF_LIBRARY})
endif()

set(separation_main_src src/separation_main.c)

if(SEPARATION_OPENCL)
  if(NVIDIA_OPENCL)
    set(CL_PLAN_CLASS "cuda_opencl")
  elseif(AMD_OPENCL)
    set(CL_PLAN_CLASS "ati_opencl")
  else()
    set(CL_PLAN_CLASS "opencl")
  endif()
else()
  set(CL_PLAN_CLASS "")
endif()

set(separation_core_libs )
# Build the critical area separately for SSE2, SSE3, no SSE2
if(SYSTEM_IS_X86 AND DOUBLEPREC AND NOT SEPARATION_CRLIBM)
  set(core_headers include/probabilities.h include/probabilities_intrin.h)
  if(HAVE_SSE2)
    add_library(separation_core_sse2 STATIC src/probabilities_sse2.c ${core_headers})
    enable_sse2(separation_core_sse2)
    disable_sse3(separation_core_sse2)
    disable_sse41(separation_core_sse2)
    list(APPEND separation_core_libs separation_core_sse2)
  endif()

  if(HAVE_SSE3)
    add_library(separation_core_sse3 STATIC src/probabilities_sse2.c ${core_headers})
    enable_sse2(separation_core_sse3)
    enable_sse3(separation_core_sse3)
    disable_sse41(separation_core_sse3)
    list(APPEND separation_core_libs separation_core_sse3)
  endif()

  if(HAVE_SSE41)
    add_library(separation_core_sse41 STATIC src/probabilities_sse2.c ${core_headers})
    enable_sse2(separation_core_sse41)
    enable_sse3(separation_core_sse41)
    enable_sse41(separation_core_sse41)
    list(APPEND separation_core_libs separation_core_sse41)
  endif()

  if(HAVE_AVX AND NOT MSVC32_AVX_WORKAROUND)
    # The real AVX one isn't done. The SSE* one is faster by building it with AVX
    #add_library(separation_core_avx STATIC src/probabilities_avx.c ${core_headers})
    add_library(separation_core_avx STATIC src/probabilities_sse2.c ${core_headers})
    enable_avx(separation_core_avx)
    list(APPEND separation_core_libs separation_core_avx)
  endif()

  if(MSVC32_AVX_WORKAROUND)
    add_definitions("-DMSVC32_AVX_WORKAROUND=1")
  endif()
endif()

add_library(separation STATIC ${lib_source_file_list} ${separation_headers})
maybe_disable_ssen(separation)


if(NOT MILKYWAY_IPHONE_APP)
  add_executable(milkyway_separation ${separation_main_src})
  maybe_disable_ssen(milkyway_separation)
else()
  add_library(milkyway_separation STATIC ${separation_main_src})
endif()

if(BOINC_RELEASE_NAMES)
  get_boinc_bin_name("milkyway_separation" "${separation_VERSION}" "${CL_PLAN_CLASS}")
  set(SEPARATION_BIN_NAME "${BOINC_BIN_NAME}" CACHE INTERNAL PARENT_SCOPE)
  set_target_properties(milkyway_separation
                          PROPERTIES
                            OUTPUT_NAME "${SEPARATION_BIN_NAME}")
endif()

milkyway_link(milkyway_separation ${BOINC_APPLICATION}
                                  ${SEPARATION_STATIC}
                                  "separation;${separation_core_libs};${exe_link_libs}")


add_subdirectory(tests EXCLUDE_FROM_ALL)

if(INSTALL_BOINC)
  install_boinc(milkyway_separation)
elseif(NOT MILKYWAY_IPHONE_APP)
  install(TARGETS milkyway_separation
            RUNTIME DESTINATION bin)
endif()

