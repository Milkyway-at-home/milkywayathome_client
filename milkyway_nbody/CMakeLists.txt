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

#No idea on the actual version, but I'll go with this for now.

cmake_minimum_required (VERSION 2.6)
project(milkyway_nbody)

include_directories("${PROJECT_SOURCE_DIR}")

set(nbody_src chisq.c
              code.c
              getparam.c
              grav.c
              io.c
              load.c
              orbitintegrator.c
              json_params.c
              util.c)

set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")

add_executable(milkyway_nbody ${nbody_src})
target_link_libraries(milkyway_nbody ${POPT_LIBRARY} ${JSON_C_LIBRARY})


