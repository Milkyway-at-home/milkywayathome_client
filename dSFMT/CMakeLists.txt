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
project(dSFMT)

set(dSFMT_VERSION_MAJOR 2 CACHE INTERNAL "dSFMT major version number")
set(dSFMT_VERSION_MINOR 1 CACHE INTERNAL "dSFMT minor version number")

set(DSFMT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}" CACHE INTERNAL "libdsfmt headers")
include_directories(${DSFMT_INCLUDE_DIR})


# Value it assumes if you don't specify it; stops warning
set(DSFMT_MEXP "19937" CACHE STRING "dSFMT mexp")
mark_as_advanced(DSFMT_MMEXP)

if(NOT MSVC)
  set(dsmft_flags "-O3 -DNDEBUG -finline-functions -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99 ${SSE2_FLAGS}")
else()
  set(dsmft_flags "${SSE2_FLAGS}")
endif()

set(dsmft_src dSFMT.c)
set(dsmft_hdr dSFMT.h)

set(DSFMT_FLAGS "-DDSFMT_MEXP=${DSFMT_MEXP}" CACHE INTERNAL "dSFMT build flags")

add_library(dsfmt STATIC ${dsmft_src})

add_definitions("${DSFMT_FLAGS}")
set_property(SOURCE ${dsmft_src}
               PROPERTY COMPILE_FLAGS "${dsmft_flags}")

# install(TARGETS dsfmt
#           ARCHIVE       DESTINATION lib
#           PUBLIC_HEADER DESTINATION include)


