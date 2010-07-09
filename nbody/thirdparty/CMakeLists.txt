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

message("arst: ${CMAKE_CURRENT_SOURCE_DIR}/dSFMT-src-${dsfmtVer}/dSFMT.c")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/dSFMT-src-${dsfmtVer}")

add_library(dsmft STATIC "${CMAKE_CURRENT_SOURCE_DIR}/dSFMT-src-${dsfmtVer}/dSFMT.c")

set(dsmft_flags "-O3 -fPIC -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -std=c99")

set_target_properties(dsmft PROPERTIES
                            COMPILE_DEFINITIONS "DSFMT_MEXP=19937"
                            COMPILE_FLAGS "${dsmft_flags} ${CMAKE_C_FLAGS}")

#TODO: The tests dSFMT has



