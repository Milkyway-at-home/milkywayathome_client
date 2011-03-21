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

find_package(OpenSSL)
include_directories(${OPENSSL_INCLUDE_DIR})

find_program(LUA_INTERPRETER lua)

add_executable(nbody_test_driver nbody_test_driver.c)
target_link_libraries(nbody_test_driver ${OPENSSL_LIBRARIES} ${nbody_exe_link_libs})


add_test(NAME checkpoint_test
         COMMAND nbody_test_driver "CheckpointTest.lua")

add_test(NAME criterion_test
         COMMAND nbody_test_driver "CriterionTest.lua")



set(body_counts 100 1024 2048)
set(seeds       670828913 886885833 715144259 430281807 543966758)

set(orphan_model_names model_1
                       model_2
                       model_3
                       model_4
                       model_5
                       model_5_bounds_test
                       model_6
                       model_7
                       model_triaxial)


foreach(n ${body_counts})
  foreach(seed ${seeds})
    foreach(model_name ${orphan_model_names})
      add_test(NAME ${model_name}__${seed}_${n}
                 COMMAND ${LUA_INTERPRETER} "RunTestUnits.lua"
                                            $<TARGET_FILE:${nbody_bin_name}>
                                            "orphan_models"
                                            ${model_name}
                                            "orphan_model_histogram"
                                            ${n}
                                            ${seed})
    endforeach()
  endforeach()
endforeach()

function(make_bodycount_test_set n)
  add_custom_target(test_${n} COMMAND ${CMAKE_CTEST_COMMAND} -R "model_.*_${n}"
                              DEPENDS ${nbody_bin_name}
                              COMMENT "Running model tests with ${n} bodies")
endfunction()

foreach(n ${body_counts})
  make_bodycount_test_set(${n})
endforeach()

add_custom_target(quick_test COMMAND ${CMAKE_CTEST_COMMAND} -R "model_.*__670828913_100"
                             DEPENDS ${nbody_bin_name})

