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

cmake_minimum_required(VERSION 2.6)
project(milkyway)

include(SpecialLibs)
include(CheckFunctionExists)
include(CheckIncludeFiles)


check_include_files(tgmath.h HAVE_TGMATH_H)
if(NOT DOUBLEPREC AND NOT HAVE_TGMATH_H)
  message(FATAL_ERROR "Float requires tgmath.h")
endif()


set(MILKYWAY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE INTERNAL "libmilkyway headers")
include_directories(${MILKYWAY_INCLUDE_DIR})

find_package(Popt REQUIRED)
include_directories(${POPT_INCLUDE_DIR})

if(BOINC_APPLICATION)
  find_package(BOINC)
  include_directories(${BOINC_INCLUDE_DIR})
endif()

if(${MILKYWAY_OPENCL} MATCHES "ON")
  find_package(OpenCL)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

set_os_specific_libs(${MILKYWAY_OPENCL})

check_function_exists(exp10 HAVE_EXP10)
check_function_exists(sincos HAVE_SINCOS)

configure_file (
  "${MILKYWAY_INCLUDE_DIR}/milkyway_config.h.in"
  "${PROJECT_INCLUDE_DIR}/milkyway_config.h")

include_directories(${MILKYWAY_INCLUDE_DIR})

set(mw_lib_src src/milkyway_util.c
               src/mw_boinc_util.c
               src/milkyway_cpp_util.cc)

set(mw_cl_lib_src src/mw_cl_show_types.c
                  src/mw_cl_device.c
                  src/mw_cl_util.c
                  src/mw_cl_program.c
                  src/mw_cl_setup.c)

set(mw_lib_headers include/milkyway_config.h.in
                   include/milkyway_extra.h
                   include/milkyway_math.h
                   include/milkyway_math_functions.h
                   include/milkyway_math_functions_cl.h
                   include/milkyway_vectors_cl.h
                   include/milkyway_math_functions_standard.h
                   include/milkyway_util.h
                   include/mw_boinc_util.h
                   include/milkyway_cpp_util.h
                   include/milkyway_vectors.h
                   include/milkyway_vectors_func.h
                   include/milkyway_vector_types.h
                   include/real.h)

set(mw_cl_lib_headers include/milkyway_cl.h
                      include/mw_cl.h
                      include/mw_cl_types.h
                      include/mw_cl_show_types.h
                      include/mw_cl_util.h
                      include/mw_cl_device.h
                      include/mw_cl_program.h
                      include/mw_cl_setup.h)

if(MILKYWAY_OPENCL)
  list(APPEND mw_lib_src ${mw_cl_lib_src})
  list(APPEND mw_lib_headers ${mw_cl_lib_headers})
endif()

add_library(milkyway STATIC ${mw_lib_src} ${mw_lib_headers})
target_link_libraries(milkyway dsfmt)


set_target_properties(milkyway PROPERTIES
                      COMPILE_DEFINITIONS "MILKYWAY_MATH_COMPILATION")


