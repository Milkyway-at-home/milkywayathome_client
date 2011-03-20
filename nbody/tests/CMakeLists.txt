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

# List of tests to run.
# Test sources should be named <testname>.c .
set(test_list t_vectors
              t_spherical
              t_miyamoto_nagai_disk
              t_exponential_disk
              t_log_halo
              t_nfw_halo
              t_triaxial_halo)


set(body_counts 100 1024 2048)
#set(seeds       670828913 886885833 715144259 430281807 543966758)
set(seeds       670828913)

set(orphan_model_names model_1
                       model_2
                       model_3
                       model_4
                       model_5
                       model_5_bounds_test
                       model_6
                       model_7
                       model_triaxial)

find_program(LUA_INTERPRETER lua)

set(nbody_test_driver "RunTestUnits.lua")
set(nbody_bin "${EXECUTABLE_OUTPUT_PATH}/${nbody_bin_name}")
set(test_histogram_name )

foreach(seed ${seeds})
  foreach(model_name ${orphan_model_names})
         add_test(${model_name}__${seed}_100 ${LUA_INTERPRETER}
                    ${nbody_test_driver} ${nbody_bin}
                                         "orphan_models"
                                         ${model_name}
                                         "orphan_model_histogram"
                                         ${seed}
                                         100)
       endforeach()
endforeach()

