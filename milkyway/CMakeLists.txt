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

check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_function_exists(asprintf HAVE_ASPRINTF)
check_include_files(tgmath.h HAVE_TGMATH_H)
check_include_files(fpu_control.h HAVE_FPU_CONTROL_H)
if(NOT DOUBLEPREC AND NOT HAVE_TGMATH_H)
  message(FATAL_ERROR "Float requires tgmath.h")
endif()


set(MILKYWAY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "libmilkyway headers")
include_directories(${MILKYWAY_INCLUDE_DIR})
include_directories(${DSFMT_INCLUDE_DIR})



if(${MILKYWAY_OPENCL} MATCHES "ON")
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

check_function_exists(exp10 HAVE_EXP10)
check_function_exists(exp2 HAVE_EXP2)
check_function_exists(sincos HAVE_SINCOS)

set(CMAKE_REQUIRED_INCLUDES "malloc.h")
check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
check_function_exists(_aligned_malloc HAVE__ALIGNED_MALLOC)
check_function_exists(__mingw_aligned_malloc HAVE___MINGW_ALIGNED_MALLOC)
set(CMAKE_REQUIRED_INCLUDES )

#check_function_exists(clock_gettime HAVE_CLOCK_GETTIME)  # WHY DOESN'T THIS WORK?
check_library_exists(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)

set(CMAKE_REQUIRED_INCLUDES "sys/time.h")
check_function_exists(gettimeofday HAVE_GETTIMEOFDAY)
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
               src/mw_boinc_util.c
               src/mw_boinc_prefs.cpp
               src/milkyway_cpp_util.cc
               src/milkyway_show.c
               src/mw_cpuid.c
               src/milkyway_timing.c)


set(milkyway_lua_src src/milkyway_lua_marshal.c
                     src/milkyway_lua_util.c
                     src/milkyway_lua_math.c
                     src/milkyway_lua_vector.c
                     src/milkyway_lua_dsfmt.c
                     src/milkyway_lua_types.c)

set(mw_cl_lib_src src/mw_cl_show_types.c
                  src/mw_cl_device.c
                  src/mw_cl_util.c
                  src/mw_cl_program.c
                  src/mw_cl_setup.c)

set(mw_lib_headers include/milkyway_config.h.in
                   include/milkyway_extra.h
                   include/milkyway_math.h
                   include/milkyway_math_functions.h
                   include/milkyway_util.h
                   include/milkyway_alloc.h
                   include/mw_boinc_util.h
                   include/mw_boinc_prefs.h
                   include/milkyway_cpp_util.h
                   include/milkyway_vectors.h
                   include/milkyway_vector_types.h
                   include/real.h
                   include/milkyway_show.h
                   include/mw_cpuid.h
                   include/milkyway_timing.h
                   include/mw_asprintf.h)


set(milkyway_lua_headers include/milkyway_lua.h
                         include/milkyway_lua_marshal.h
                         include/milkyway_lua_util.h
                         include/milkyway_lua_math.h
                         include/milkyway_lua_vector.h
                         include/milkyway_lua_dsfmt.h
                         include/milkyway_lua_types.h)

set(mw_cl_lib_headers include/milkyway_cl.h
                      include/mw_cl.h
                      include/mw_cl_types.h
                      include/mw_cl_show_types.h
                      include/mw_cl_util.h
                      include/mw_cl_device.h
                      include/mw_cl_program.h
                      include/mw_cl_setup.h)

if(WIN32 AND NOT HAVE_ASPRINTF)
    list(APPEND mw_lib_src src/mw_asprintf.c)
endif()



add_library(milkyway_lua STATIC ${milkyway_lua_src} ${milkyway_lua_headers})
set_target_properties(milkyway_lua PROPERTIES
                        COMPILE_DEFINITIONS "MILKYWAY_LUA_COMPILATION")

target_link_libraries(milkyway_lua milkyway dsfmt)

add_library(milkyway STATIC ${mw_lib_src} ${mw_lib_headers})
target_link_libraries(milkyway dsfmt)
set_target_properties(milkyway PROPERTIES
                        COMPILE_DEFINITIONS "MILKYWAY_MATH_COMPILATION")

# install(TARGETS milkyway milkyway_lua
#         ARCHIVE       DESTINATION lib
#         PUBLIC_HEADER DESTINATION include)


if(MILKYWAY_OPENCL)
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

