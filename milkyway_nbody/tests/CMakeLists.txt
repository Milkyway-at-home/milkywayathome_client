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


include_directories("${PROJECT_SOURCE_DIR}")

# List of tests to run.
# Test sources should be named <testname>.c .
set(test_list t_vectors
              t_spherical
              t_loghalo
              t_miyamoto_nagai
              )

set(nbody_main_src ../plummer.c
                   ../nbody_simulation.c
                   ../grav.c
                   ../load.c
                   ../util.c
                   ../io.c
                   ../json_params.c
                   ../orbitintegrator.c)


foreach(test ${test_list})
	add_executable(${test} ${test}.c ${nbody_main_src})
	add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
    add_dependencies(check ${test})
    target_link_libraries(${test} m ${JSON_C_LIBRARY})
endforeach()

