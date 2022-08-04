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
project(milkyway)

include(SpecialLibs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(MWLink)
include(CheckSSE2)
include(GetGitRevisionDescription)

get_git_head_revision(GIT_REFSPEC MILKYWAY_GIT_COMMIT_ID)
git_describe(MILKYWAY_GIT_DESCRIBE)

check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(windows.h HAVE_WINDOWS_H)
check_include_files(direct.h HAVE_DIRECT_H)
check_include_files(fcntl.h HAVE_FCNTL_H)
check_include_files(signal.h HAVE_SIGNAL_H)
check_include_files(process.h HAVE_PROCESS_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(malloc.h HAVE_MALLOC_H)
check_function_exists(asprintf HAVE_ASPRINTF)
check_include_files(float.h HAVE_FLOAT_H)
check_include_files(fpu_control.h HAVE_FPU_CONTROL_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(sys/mman.h HAVE_SYS_MMAN_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/wait.h HAVE_SYS_WAIT_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)

set(MILKYWAY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "libmilkyway headers")
include_directories(${MILKYWAY_INCLUDE_DIR})
include_directories(${DSFMT_INCLUDE_DIR})



if(${MILKYWAY_OPENCL} MATCHES "ON")
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(exp10 HAVE_EXP10)
check_function_exists(exp2 HAVE_EXP2)
check_function_exists(sincos HAVE_SINCOS)
check_function_exists(fmax HAVE_FMAX)
check_function_exists(fmin HAVE_FMIN)
set(CMAKE_REQUIRED_LIBRARIES )


set(CMAKE_REQUIRED_INCLUDES "malloc.h")
check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
check_function_exists(_aligned_malloc HAVE__ALIGNED_MALLOC)
check_function_exists(__mingw_aligned_malloc HAVE___MINGW_ALIGNED_MALLOC)
set(CMAKE_REQUIRED_INCLUDES )

set(CMAKE_REQUIRED_LIBRARIES rt)
check_function_exists(clock_gettime HAVE_CLOCK_GETTIME)
set(CMAKE_REQUIRED_LIBRARIES )

set(CMAKE_REQUIRED_INCLUDES "sys/time.h")
check_function_exists(gettimeofday HAVE_GETTIMEOFDAY)
set(CMAKE_REQUIRED_INCLUDES )

set(CMAKE_REQUIRED_INCLUDES "signal.h")
check_function_exists(sigaction HAVE_SIGACTION)
set(CMAKE_REQUIRED_INCLUDES )


set(CMAKE_REQUIRED_INCLUDES "mach/mach_time.h")
check_function_exists(mach_absolute_time HAVE_MACH_ABSOLUTE_TIME)
set(CMAKE_REQUIRED_INCLUDES )



configure_file(
  "${MILKYWAY_INCLUDE_DIR}/milkyway_config.h.in"
  "${PROJECT_INCLUDE_DIR}/milkyway_config.h")

configure_file(
  "${MILKYWAY_INCLUDE_DIR}/milkyway_git_version.h.in"
  "${PROJECT_INCLUDE_DIR}/milkyway_git_version.h")


include_directories(${MILKYWAY_INCLUDE_DIR})

# Careful with adding math functions here Right now the crlibm setting
# is individual for separation vs. nbody. Math here won't be compiled
# with crlibm, leading to unexpected results
set(mw_lib_src src/milkyway_util.c
               src/milkyway_alloc.c
               src/milkyway_rename.c
               src/milkyway_boinc_util.cc
               src/milkyway_show.c
               src/milkyway_cpuid.c
               src/milkyway_timing.c)


set(milkyway_lua_src src/milkyway_lua_marshal.c
                     src/milkyway_lua_util.c
                     src/milkyway_lua_math.c
                     src/milkyway_lua_vector.c
                     src/milkyway_lua_dsfmt.c
                     src/milkyway_lua_types.c)

set(mw_cl_lib_src src/milkyway_cl_show_types.c
                  src/milkyway_cl_device.c
                  src/milkyway_cl_util.c
                  src/milkyway_cl_program.c
                  src/milkyway_cl_setup.c)

set(mw_lib_headers include/milkyway_config.h.in
                   include/milkyway_extra.h
                   include/milkyway_math.h
                   include/milkyway_math_double.h
                   include/milkyway_math_float.h
                   include/milkyway_math_supplemental.h
                   include/milkyway_util.h
                   include/milkyway_alloc.h
                   include/milkyway_boinc_util.h
                   include/milkyway_rename.h
                   include/milkyway_vectors.h
                   include/milkyway_show.h
                   include/milkyway_cpuid.h
                   include/milkyway_timing.h
                   include/milkyway_asprintf.h
                   include/milkyway_simd_defs.h
                   include/milkyway_sse2_intrin.h
                   include/milkyway_simd.h)


set(milkyway_lua_headers include/milkyway_lua.h
                         include/milkyway_lua_marshal.h
                         include/milkyway_lua_util.h
                         include/milkyway_lua_math.h
                         include/milkyway_lua_vector.h
                         include/milkyway_lua_dsfmt.h
                         include/milkyway_lua_types.h)

set(mw_cl_lib_headers include/milkyway_cl.h
                      include/mw_cl.h
                      include/milkyway_cl_types.h
                      include/milkyway_cl_show_types.h
                      include/milkyway_cl_util.h
                      include/milkyway_cl_device.h
                      include/milkyway_cl_program.h
                      include/milkyway_cl_setup.h)

if(WIN32 AND NOT HAVE_ASPRINTF)
    list(APPEND mw_lib_src src/milkyway_asprintf.c)
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR}/../boinc/boinc)

add_library(milkyway_lua STATIC ${milkyway_lua_src} ${milkyway_lua_headers})
set_target_properties(milkyway_lua PROPERTIES
                        COMPILE_DEFINITIONS "MILKYWAY_LUA_COMPILATION")

target_link_libraries(milkyway_lua milkyway dsfmt)

add_library(milkyway STATIC ${mw_lib_src} ${mw_lib_headers})
target_link_libraries(milkyway dsfmt)
if(BOINC_APPLICATION)
  target_link_libraries(milkyway ${BOINC_LIBRARIES})
endif()
set_target_properties(milkyway PROPERTIES
                        COMPILE_DEFINITIONS "MILKYWAY_MATH_COMPILATION")

# install(TARGETS milkyway milkyway_lua
#         ARCHIVE       DESTINATION lib
#         PUBLIC_HEADER DESTINATION include)


if(MILKYWAY_OPENCL)
  add_definitions("-DCL_USE_DEPRECATED_OPENCL_1_1_APIS=1")   # For now stick to 1.1
  add_library(milkyway_cl STATIC ${mw_cl_lib_src} ${mw_cl_lib_headers})
  target_link_libraries(milkyway_cl milkyway)
  set_target_properties(milkyway_cl PROPERTIES
                          COMPILE_DEFINITIONS "MILKYWAY_CL_COMPILATION")


# install(TARGETS milkyway_cl
#         ARCHIVE       DESTINATION lib
#         PUBLIC_HEADER DESTINATION include)
endif()



if(WIN32)
  target_link_libraries(milkyway winmm)
endif()

if(HAVE_LIBRT)
  target_link_libraries(milkyway rt)
endif()

maybe_disable_ssen(milkyway milkyway_lua)
if(MILKYWAY_OPENCL)
    maybe_disable_ssen(milkyway_cl)
endif()

