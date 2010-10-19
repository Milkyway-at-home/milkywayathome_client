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

cmake_minimum_required (VERSION 2.6)
project(milkyway_nbody)

set(nbody_MAJOR_VERSION 0)
set(nbody_MINOR_VERSION 21)
set(nbody_VERSION "${nbody_MAJOR_VERSION}.${nbody_MINOR_VERSION}")


include(MaybeDlCheck)
include(StaticBOINC)
include(SpecialLibs)
include(PrintInfo)
include(CheckFunctionExists)

maybe_static()

if(BOINC_APPLICATION)
  find_package(BOINC)
  include_directories(${BOINC_INCLUDE_DIR})
endif()

find_package(JsonC REQUIRED)
include_directories(${JSON_C_INCLUDE_DIR})

check_function_exists(asprintf HAS_ASPRINTF)

if(NEED_POPT_DEPS)
    #These are only dependencies of popt needed for static linking
  find_package(Libintl REQUIRED)
  find_package(Libiconv REQUIRED)
endif()
find_package(Popt REQUIRED)
include_directories(${POPT_INCLUDE_DIR})

set(NBODY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
include_directories(${NBODY_INCLUDE_DIR})

option(BOINC_APPLICATION
  "Build as a BOINC Application" ON)

option(BOINC_DEBUG
  "Avoid annoying things BOINC does for debugging" OFF)

option(NBODY_OPENCL
       "Use OpenCL" OFF)

option(ENABLE_CRLIBM
      "Use crlibm double" ON)
mark_as_advanced(ENABLE_CRLIBM)

option(NBODY_OPENMP "Use OpenMP" OFF)

find_package(OpenMP)

if(OPENMP_FOUND AND NBODY_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
elseif(NOT OPENMP_POUND AND NBODY_OPENMP)
  message(WARNING "Did not find OpenMP support, but was enabled. Continuing without OpenMP")
endif()



if(BOINC_RELEASE)
  set(ENABLE_CRLIBM ON)
endif()

if(NOT DOUBLEPREC AND ENABLE_CRLIBM)
  message(WARNING "Attempting to use single precision with crlibm. Disabling use of crlibm")
  set(ENABLE_CRLIBM OFF)
endif()

if(NBODY_OPENCL)
  find_package(OpenCL)
  if(NOT OPENCL_FOUND)
    message(WARNING "N-body OpenCL was requested, but OpenCL not found. Disabling OpenCL.")
    set(NBODY_OPENCL OFF)
  else()
    include_directories(${OPENCL_DIRECTORIES})
  endif()
endif()

set_os_specific_libs(${NBODY_OPENCL})


# Handle downloading, extracting and building dSFMT

set(dsfmtVer "2.1")
set(randSrcMD5 "b3a38dac7fd8996a70d02edc4432dd75")
set(randSrcPath "${MILKYWAY_THIRDPARTY}")
set(dsfmtTar "dSFMT-src-${dsfmtVer}.tar.gz")
maybe_dl_check("dSFMT"
               "${dsfmtVer}"
               "${randSrcMD5}"
               "http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/${dsfmtTar}"
               "${randSrcPath}"
               "${dsfmtTar}")


set(dsfmt_mexp "19937") #Value it assumes if you don't specify it; stops warning
set(dsmft_flags "-O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99 -mfpmath=sse -msse -msse2")

set(dsmft_src "${MILKYWAY_THIRDPARTY}/dSFMT-src-${dsfmtVer}/dSFMT.c")
include_directories("${MILKYWAY_THIRDPARTY}/dSFMT-src-${dsfmtVer}/")
add_definitions("-DDSFMT_MEXP=${dsfmt_mexp}")



print_separator()
message("   Building:          ${PROJECT_NAME}")
message("   BOINC Application: ${BOINC_APPLICATION}")
message("   BOINC Debug:       ${BOINC_DEBUG}")
message("   Double precision:  ${DOUBLEPREC}")
message("   Enable crlibm:     ${ENABLE_CRLIBM}")
message("  ----")
print_libs()
print_separator()

if(BOINC_APPLICATION AND NOT BOINC_LIBRARIES)
  message(FATAL_ERROR "Trying to build as a BOINC application, but the BOINC library wasn't found")
endif()



#CHECKME: We don't actually need special number math
#except for assignment in a couple places. Not sure how useful this is
#anyway.
#if(NOT ENABLE_FDLIBM)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffinite-math-only")
#endif()

set(NBODY_SRC_DIR "${PROJECT_SOURCE_DIR}/src/")
set(NBODY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include/")
set(nbody_lib_src ${NBODY_SRC_DIR}/chisq.c
                  ${NBODY_SRC_DIR}/grav.c
                  ${NBODY_SRC_DIR}/io.c
                  ${NBODY_SRC_DIR}/load.c
                  ${NBODY_SRC_DIR}/orbitintegrator.c
                  ${NBODY_SRC_DIR}/accelerations.c
                  ${NBODY_SRC_DIR}/nbody_simulation.c
                  ${NBODY_SRC_DIR}/json_params.c
                  ${NBODY_SRC_DIR}/calc_params.c
                  ${NBODY_SRC_DIR}/plummer.c
                  ${NBODY_SRC_DIR}/nbody_util.c
                  ${NBODY_SRC_DIR}/nbody_step.c
                  ${dsmft_src})

set(nbody_lib_headers ${NBODY_INCLUDE_DIR}/chisq.h
                      ${NBODY_INCLUDE_DIR}/grav.h
                      ${NBODY_INCLUDE_DIR}/nbody_config.h.in
                      ${NBODY_INCLUDE_DIR}/io.h
                      ${NBODY_INCLUDE_DIR}/load.h
                      ${NBODY_INCLUDE_DIR}/orbitintegrator.h
                      ${NBODY_INCLUDE_DIR}/accelerations.h
                      ${NBODY_INCLUDE_DIR}/json_params.h
                      ${NBODY_INCLUDE_DIR}/calc_params.h
                      ${NBODY_INCLUDE_DIR}/plummer.h
                      ${NBODY_INCLUDE_DIR}/nbody.h
                      ${NBODY_INCLUDE_DIR}/nbody_priv.h
                      ${NBODY_INCLUDE_DIR}/nbody_types.h
                      ${NBODY_INCLUDE_DIR}/nbody_util.h)

if(BOINC_RELEASE OR NOT HAS_ASPRINTF)
  set(DEFAULT_NBODY_SHOW OFF)
else()
  set(DEFAULT_NBODY_SHOW ON)
endif()

option(NBODY_SHOW "Build show type functions" ${DEFAULT_NBODY_SHOW})

if(NBODY_SHOW)
  #Uses asprintf extension
  list(APPEND nbody_lib_src ${NBODY_SRC_DIR}/show.c)
  list(APPEND nbody_lib_headers ${NBODY_INCLUDE_DIR}/show.h)
endif()


if(BOINC_APPLICATION)
  list(APPEND nbody_lib_src ${NBODY_SRC_DIR}/checkpoint.c)
  list(APPEND nbody_lib_headers ${NBODY_INCLUDE_DIR}/checkpoint.h)
endif()

set_property(SOURCE ${dsmft_src}
               PROPERTY COMPILE_DEFINITIONS "DSFMT_MEXP=19937")

set_property(SOURCE ${dsmft_src}
               PROPERTY COMPILE_FLAGS "${dsmft_flags}")

#TODO: The tests dSFMT has



# Require SSE and SSE2 (needed for doubles) This avoids a host of
# issues with x87, ranging from different intermediate precisions on
# different places, and removes many of the places where different
# compiler optimizations can change the precision.
#
# Some of these issues are summarized here: There are others to be
# aware of.
# http://www.yosefk.com/blog/consistency-how-to-defeat-the-purpose-of-ieee-floating-point.html
#
# Also as a nice advantage, SSE for a quick benchmark I ran is ~20%
# faster than using the antique x87 junk. There is also a later
# opportunity to do hand vectorization for an even bigger improvement.
#
# The GPU floating point stuff is
# "similar" to SSE, but not SSE.  OpenCL demands a stricter
# conformance to standards for doubles, so hopefully when we start
# using GPUs the results will be the same as using SSE/SSE2 and
# doubles. No promises however.
#
# Summary: we abandon old, old processors for the sake of convenience
# and consistency, and a nice speed boost on modern x86/x86_64
# hardware.
#

set(CMAKE_C_FLAGS "-std=gnu99 -mfpmath=sse -msse -msse2 ${CMAKE_C_FLAGS}")

configure_file (
  "${NBODY_INCLUDE_DIR}/nbody_config.h.in"
  "${PROJECT_INCLUDE_DIR}/nbody_config.h")

if(WIN32)
  add_definitions(-D_WIN32_IE=0x0501)
  add_definitions(-D_WINDOWS)
  add_definitions(-D_WIN32_WINDOWS=0x0501)
endif()

list(APPEND nbody_link_libs ${JSON_C_LIBRARY} ${OS_SPECIFIC_LIBS} milkyway ${LIBM_LIBRARY})

if(NBODY_OPENCL)
  set(c_kernel_dir "${CMAKE_CURRENT_SOURCE_DIR}/ckernels")
  list(APPEND nbody_lib_src "${PROJECT_SOURCE_DIR}/src/gravmap_opencl.c")

  list(APPEND nbody_lib_src "${c_kernel_dir}/cl_gravmap.c")
  set_source_files_properties("${c_kernel_dir}/cl_gravmap.c"
                              PROPERTIES GENERATED TRUE)

  list(APPEND nbody_link_libs ${OPENCL_LIBRARIES})
endif()


add_library(nbody_simulation STATIC ${nbody_lib_src} ${nbody_lib_headers})

set(nbody_libs nbody_simulation)

set(exe_link_libs ${LIBM_LIBRARY}
                  ${POPT_LIBRARY}
                  ${nbody_libs}
                  ${nbody_link_libs})


if(ENABLE_CRLIBM)
  if(DOWNLOAD_CRLIBM)
    add_dependencies(nbody_simulation crlibm_build)
    set(CRLIBM_INCLUDE_DIR "${MILKYWAY_ROOT}/include"})
    set(CRLIBM_LIBRARY "${LIBRARY_OUTPUT_PATH}/libcrlibm${CMAKE_STATIC_LIBRARY_SUFFIX}")
  else()
    if(BOINC_RELEASE)
      find_package(Crlibm REQUIRED)
    else()
      find_package(Crlibm)
    endif()
  endif()

  include_directories(${CRLIBM_INCLUDE_DIR})
  list(APPEND exe_link_libs ${CRLIBM_LIBRARY})
endif()

find_package(Libm REQUIRED)

if(NEED_POPT_DEPS)
  # Popt dependencies for static linking For some reason, on OS X /
  # windows you need these dependencies but not Linux.
  list(APPEND exe_link_libs ${LIBICONV_LIBRARY}
                            ${LIBINTL_LIBRARY})
endif()



if(BOINC_RELEASE_NAME)
  if(NBODY_OPENMP)
    set(NBODY_PLAN "mt")
  else()
    set(NBODY_PLAN "sse2")
  endif()

  get_boinc_bin_name("milkyway_nbody" "${nbody_VERSION}" "${NBODY_PLAN}")
  set(nbody_bin_name "${BOINC_BIN_NAME}")
else()
  set(nbody_bin_name "milkyway_nbody")
endif()



add_executable(${nbody_bin_name} ${NBODY_SRC_DIR}/main.c)
milkyway_link(${nbody_bin_name} ${BOINC_APPLICATION} ${BOINC_RELEASE} "${exe_link_libs}")
#add_subdirectory(tests EXCLUDE_FROM_ALL)


