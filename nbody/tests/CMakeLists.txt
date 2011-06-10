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

include_directories(${OPENSSL_INCLUDE_DIR})

add_executable(nbody_test_driver nbody_test_driver.c)
target_link_libraries(nbody_test_driver ${OPENSSL_LIBRARIES} ${nbody_exe_link_libs})

if(BOINC_APPLICATION)
  if(UNIX)
    target_link_libraries(nbody_test_driver pthread)
  endif()
endif()


add_test(NAME checkpoint_test
         COMMAND nbody_test_driver "CheckpointTest.lua")

# add_test(NAME context_test
#          COMMAND nbody_test_driver "ContextTest.lua")

# Run these first 3 tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -I 1,2,1)

# TODO: Test with known tree incest to check return status right


set(body_counts 100 1024 10000)

set(orphan_model_names model_1
                       model_2
                       model_3
                       model_4
                       model_5
                       model_5_bounds_test
                       model_6
                       model_7
                       model_triaxial)

foreach(model_name ${orphan_model_names})
  foreach(n ${body_counts})
    add_test(NAME ${model_name}__${n}_test
               COMMAND nbody_test_driver "RunTestUnits.lua"
                                          $<TARGET_FILE:milkyway_nbody>
                                          "orphan_models"
                                          ${model_name}
                                          "orphan_model_histogram"
                                          ${n})
  endforeach()
endforeach()

function(make_bodycount_test_set n)
  add_custom_target(test_${n} COMMAND ${CMAKE_CTEST_COMMAND} -R "model_.*__${n}_test$"
                              DEPENDS ${nbody_bin_name}
                              COMMENT "Running model tests with ${n} bodies")
endfunction()

foreach(n ${body_counts})
  make_bodycount_test_set(${n})
endforeach()


# compress: ${CMAKE_COMMAND} -E tar cvfj test_output.tar.bz2 test_output
#execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvfj "context_test_results.tar.bz2"
#                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")


add_custom_target(context_test COMMAND ${CMAKE_CTEST_COMMAND} -R "context_test")



