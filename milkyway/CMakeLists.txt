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

set(MILKYWAY_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
include_directories(${MILKYWAY_INCLUDE_DIR})

if(BOINC_APPLICATION)
  find_package(BOINC)
  include_directories(${BOINC_INCLUDE_DIR})
endif()


configure_file (
  "${MILKYWAY_INCLUDE_DIR}/milkyway_config.h.in"
  "${PROJECT_INCLUDE_DIR}/milkyway_config.h")

set(mw_lib_src src/milkyway_util.c)
set(mw_cl_lib_src src/show_cl_types.c
                  src/build_cl.c)

if(MILKYWAY_OPENCL)
  list(APPEND mw_lib_src ${mw_cl_lib_src})
endif()

add_library(milkyway STATIC ${mw_lib_src})

set_target_properties(milkyway PROPERTIES
                      COMPILE_DEFINITIONS "MILKYWAY_MATH_COMPILATION")


