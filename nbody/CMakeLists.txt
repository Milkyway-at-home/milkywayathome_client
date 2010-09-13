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

cmake_minimum_required (VERSION 2.8)
project(milkyway_nbody)

set(milkyway_nbody_MAJOR_VERSION 0)
set(milkyway_nbody_MINOR_VERSION 04)
set(milkyway_nbody_VERSION "${milkyway_nbody_MAJOR_VERSION}.${milkyway_nbody_MINOR_VERSION}")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(CMAKE_SKIP_BUILD_RPATH TRUE)

include(MaybeDlCheck)
include(StaticBOINC)

include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${LIBRARY_OUTPUT_PATH})

option(BOINC_APPLICATION
  "Build as a BOINC Application" ON)

option(BOINC_DEBUG
  "Avoid annoying things BOINC does for debugging" OFF)

#option(NBODY_OPENCL
#  "Use OpenCL" ${OPENCL_FOUND})

option(NBODY_OPENCL
     "Use OpenCL" OFF)


option(ENABLE_CRLIBM
      "Use crlibm double" ON)
mark_as_advanced(ENABLE_CRLIBM)

if(NOT DOUBLEPREC)
  set(ENABLE_CRLIBM OFF)
endif()


set(nbody_thirdparty "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty")

#This is really only intended for windows, since it's a pain in the
#ass. Everywhere else has much more convenient ways of installing.
if(WIN32)
  set(default_use_system_json_c OFF)
else()
  set(default_use_system_json_c ON)
endif()

option(USE_SYSTEM_JSON_C
       "Use system json-c" ${default_use_system_json_c})

if(NOT USE_SYSTEM_JSON_C)

  find_path(JSON_C_INCLUDE_DIR "json/json.h" "${PROJECT_INCLUDE_DIR}")
  find_library(JSON_C_LIBRARY json "${LIBRARY_OUTPUT_PATH}")

  # set(jsoncVer "0.9")
  # set(jsoncMD5 "3a13d264528dcbaf3931b0cede24abae")
  # set(jsoncSrcPath "${nbody_thirdparty}")
  # set(jsoncTar "json-c-${jsoncVer}.tar.gz")
  # maybe_dl_check("json-c"
  #                "${jsoncVer}"
  #                "${jsoncMD5}"
  #                "http://oss.metaparadigm.com/json-c/${jsoncTar}"
  #                "${jsoncSrcPath}"
  #                "${jsoncTar}")

  # set(jsonc_file "${LIBRARY_OUTPUT_PATH}/libjson${CMAKE_STATIC_LIBRARY_SUFFIX}")

  # if(WIN32)
  #   file(RENAME "${jsoncSrcPath}/json-c-${crlibmVer}/config.h.win32"
  #               "${jsoncSrcPath}/json-c-${crlibmVer}/config.h")
  # endif()

  # add_custom_command(
  #   OUTPUT "${jsonc_file}"
  #   COMMAND "${jsoncSrcPath}/json-c-${jsoncVer}/configure"
  #                        "--prefix=${CMAKE_CURRENT_SOURCE_DIR}"
  #                        "--disable-werror"
  #                        "CC=${CMAKE_C_COMPILER}"
  #                        "MAKE=${CMAKE_BUILD_TOOL}"
  #                        "CFLAGS=${CMAKE_C_FLAGS}"
  #   COMMAND "${CMAKE_BUILD_TOOL}"
  #   COMMAND "${CMAKE_BUILD_TOOL}" "install"
  #   WORKING_DIRECTORY "${jsoncSrcPath}/json-c-${jsoncVer}/"
  #   COMMENT "Building json-c")

  # set(JSON_C_LIBRARY ${jsonc_file})

  # add_custom_target(jsonc_build DEPENDS "${jsonc_file}")
else()
  find_package(JsonC)
  include_directories(${JSON_C_INCLUDE_DIR})
endif()


# Handle downloading, extracting and building dSFMT

set(dsfmtVer "2.1")
set(randSrcMD5 "b3a38dac7fd8996a70d02edc4432dd75")
set(randSrcPath "${nbody_thirdparty}")
set(dsfmtTar "dSFMT-src-${dsfmtVer}.tar.gz")
maybe_dl_check("dSFMT"
               "${dsfmtVer}"
               "${randSrcMD5}"
               "http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/${dsfmtTar}"
               "${randSrcPath}"
               "${dsfmtTar}")


set(dsfmt_mexp "19937") #Value it assumes if you don't specify it; stops warning
set(dsmft_flags "-O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99 -mfpmath=sse -msse -msse2")

#Find the headers for dSFMT


include_directories("${nbody_thirdparty}/dSFMT-src-${dsfmtVer}/")

set(dsmft_src "${nbody_thirdparty}/dSFMT-src-${dsfmtVer}/dSFMT.c")
add_definitions("-DDSFMT_MEXP=${dsfmt_mexp}")


# Take care of crlibm
if(ENABLE_CRLIBM)
  set(crlibmVer "1.0beta4")
  set(crlibmMD5 "8ecabd55d4a2d08030eb770ae2a5670a")
  set(crlibmTar "crlibm-${crlibmVer}.tar.gz")
  set(crlibmSrcPath "${nbody_thirdparty}")

  maybe_dl_check("crlibm"
                 "${crlibmVer}"
                 "${crlibmMD5}"
                 "http://lipforge.ens-lyon.fr/frs/download.php/162/${crlibmTar}"
                 "${nbody_thirdparty}"
                 "${crlibmTar}")
  #add_custom_target(crlibm ALL)

  #TODO: Set march, mtune etc. from user flags, but avoid some others
  set(crlibm_flags "-std=c99 -O3 -mfpmath=sse -msse -msse2")

  set(crlibm_file "${LIBRARY_OUTPUT_PATH}/libcrlibm${CMAKE_STATIC_LIBRARY_SUFFIX}")
  set(crlibm_header "${PROJECT_INCLUDE_DIR}/crlibm.h")
  set(scs_file "${LIBRARY_OUTPUT_PATH}/libscs{CMAKE_STATIC_LIBRARY_SUFFIX}")
  add_custom_command(
        OUTPUT "${crlibm_file}" "${crlibm_header}"
        COMMAND "${crlibmSrcPath}/crlibm-${crlibmVer}/configure"
                                                 "--enable-sse2"
                                                 "--prefix=${CMAKE_CURRENT_SOURCE_DIR}"
                                                 "CC=${CMAKE_C_COMPILER}"
                                                 "MAKE=${CMAKE_BUILD_TOOL}"
                                                 "CFLAGS=${crlibm_flags}"
        COMMAND "${CMAKE_BUILD_TOOL}"
        COMMAND "${CMAKE_BUILD_TOOL}" "install"
        WORKING_DIRECTORY "${crlibmSrcPath}/crlibm-${crlibmVer}"
        COMMENT "Building crlibm")

  add_custom_target(crlibm_build DEPENDS "${crlibm_file}")
  set(nbody_link_libs crlibm)
endif()


#CHECKME: Seems to only be problematic when we define our own NAN, at
#least for GCC, since we don't actually need special number math
#except for assignment in a couple places. Not sure how useful this is
#anyway.
#if(NOT ENABLE_FDLIBM)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffinite-math-only")
#endif()

set(nbody_lib_src ${PROJECT_SOURCE_DIR}/src/chisq.c
                  ${PROJECT_SOURCE_DIR}/src/grav.c
                  ${PROJECT_SOURCE_DIR}/src/io.c
                  ${PROJECT_SOURCE_DIR}/src/load.c
                  ${PROJECT_SOURCE_DIR}/src/orbitintegrator.c
                  ${PROJECT_SOURCE_DIR}/src/accelerations.c
                  ${PROJECT_SOURCE_DIR}/src/nbody_simulation.c
                  ${PROJECT_SOURCE_DIR}/src/json_params.c
                  ${PROJECT_SOURCE_DIR}/src/show.c
                  ${PROJECT_SOURCE_DIR}/src/plummer.c
                  ${PROJECT_SOURCE_DIR}/src/nbody_util.c
                  ${dsmft_src})

if(BOINC_APPLICATION)
  list(APPEND nbody_lib_src ${PROJECT_SOURCE_DIR}/src/checkpoint.c)
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

if(APPLE)
  # Try to avoid the dyld: unknown required load command 0x80000022
  # runtime error on Leopard for binaries built on 10.6
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isysroot/Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5")
endif()


message("---------------------------------------------")
message("   Building:          ${PROJECT_NAME}")
message("   BOINC Application: ${BOINC_APPLICATION}")
message("   BOINC Debug:       ${BOINC_DEBUG}")
message("   Double precision:  ${DOUBLEPREC}")
message("   Enable crlibm:     ${ENABLE_CRLIBM}")
message("   json-c Library     ${JSON_C_LIBRARY}")
message("   BOINC Libraries    ${BOINC_LIBRARIES}")
message("   Popt Library       ${POPT_LIBRARY}")
message("---------------------------------------------")

if(BOINC_APPLICATION AND NOT BOINC_LIBRARIES)
  message(FATAL_ERROR "Trying to build as a BOINC application, but the BOINC library wasn't found")
endif()


configure_file (
  "${PROJECT_INCLUDE_DIR}/nbody_config.h.in"
  "${PROJECT_INCLUDE_DIR}/nbody_config.h")

if(WIN32)
  list(APPEND OS_SPECIFIC_LIBS "msvcrt" "stdc++")
  add_definitions(-D_WIN32_IE=0x0501)
  add_definitions(-D_WINDOWS)
  add_definitions(-D_WIN32_WINDOWS=0x0501)
endif()

list(APPEND nbody_link_libs ${JSON_C_LIBRARY} ${OS_SPECIFIC_LIBS} m)

if(BOINC_APPLICATION)
  list(APPEND nbody_link_libs ${BOINC_LIBRARIES})
endif()

if(NBODY_OPENCL)
  set(c_kernel_dir "${CMAKE_CURRENT_SOURCE_DIR}/ckernels")
  list(APPEND nbody_lib_src "${PROJECT_SOURCE_DIR}/src/gravmap_opencl.c")

  list(APPEND nbody_lib_src "${c_kernel_dir}/cl_gravmap.c")
  set_source_files_properties("${c_kernel_dir}/cl_gravmap.c"
                              PROPERTIES GENERATED TRUE)

  list(APPEND nbody_link_libs ${OPENCL_LIBRARIES})
endif()


add_library(nbody_simulation STATIC ${nbody_lib_src})

set(nbody_libs nbody_simulation)

#Inline the CL kernels if they have changed, and make the library
#depend on them
if(NBODY_OPENCL)
  add_custom_command(OUTPUT "${c_kernel_dir}/cl_gravmap.c" "${c_kernel_dir}/cl_gravmap.h"
                     COMMAND ${CMAKE_COMMAND}
                               "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}"
                               "-DKERNEL_NAME=gravmap"
                               "-DKERNEL_FILE=${CMAKE_CURRENT_SOURCE_DIR}/kernels/gravmap.cl"
                               "-DC_KERNEL_DIR=${CMAKE_CURRENT_SOURCE_DIR}/ckernels"
                               -P "${CMAKE_MODULE_PATH}/InlineNBody.cmake"
                      COMMENT "Inlining gravmap kernel")

  # I have no idea why you need to do this extra target. I'm probably
  # doing it wrong.
  add_custom_target(gravmap_kernel DEPENDS "${c_kernel_dir}/cl_gravmap.c")

  add_dependencies(nbody_simulation gravmap_kernel)
endif()

if (ENABLE_CRLIBM)
  add_dependencies(nbody_simulation crlibm_build)
endif()

if(NOT USE_SYSTEM_JSON_C)
  add_dependencies(nbody_simulation jsonc_build)
endif()

add_executable(milkyway_nbody ${PROJECT_SOURCE_DIR}/src/main.c)

set(exe_link_libs ${POPT_LIBRARY}
                  ${nbody_libs}
                  ${nbody_link_libs}
                  milkyway)

if(NEED_POPT_DEPS)
  # Popt dependencies for static linking For some reason, on OS X /
  # windows you need these dependencies but not Linux.
  list(APPEND exe_link_libs ${LIBICONV_LIBRARY}
                            ${LIBINTL_LIBRARY})
endif()

correct_static_link(milkyway_nbody)
if(BOINC_APPLICATION)
  list(APPEND exe_link_libs ${BOINC_LIBRARIES})
endif()

target_link_libraries(milkyway_nbody ${exe_link_libs})

#add_subdirectory(tests EXCLUDE_FROM_ALL)


