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


cmake_dependent_option(USE_SSL_TESTS "Use tests that use OpenSSL" ON
                                      "OPENSSL_FOUND" OFF)
mark_as_advanced(USE_SSL_TESTS)

if(USE_SSL_TESTS)
  include_directories(${OPENSSL_INCLUDE_DIR})
  add_definitions("-DUSE_SSL_TESTS=1")
  set(test_driver_libs "${OPENSSL_LIBRARIES};${nbody_exe_link_libs}")
else()
  add_definitions("-DUSE_SSL_TESTS=0")
  set(test_driver_libs "${nbody_exe_link_libs}")
endif()


add_executable(nbody_test_driver nbody_test_driver.c)
milkyway_link(nbody_test_driver ${BOINC_APPLICATION} ${NBODY_STATIC} "${test_driver_libs}")

add_executable(emd_test emd_test.c)

set(emd_test_link_libs nbody
                       milkyway)

add_executable(bessel_test bessel_test.c)

set(bessel_test_link_libs nbody
                          milkyway)

add_executable(poisson_test poisson_test.c)

set(poisson_test_link_libs "${nbody_exe_link_libs}")

add_executable(virial_test virial_test.c)

set(virial_test_link_libs nbody
                          milkyway)

add_executable(average_bins_test average_bins_test.c)

set(average_bins_test_link_libs "${nbody_exe_link_libs}")

add_executable(propermotion_test propermotion_test.c)

set(propermotion_test_link_libs "${nbody_exe_link_libs}")

add_executable(EMD_Range_test EMD_Range_test.c)

set(EMD_Range_test_link_libs "${nbody_exe_link_libs}")

add_executable(mixeddwarf_test mixeddwarf_test.c test_env_util.c)

set(mixeddwarf_test_link_libs "${nbody_exe_link_libs}")

if(NBODY_CRLIBM)
    list(APPEND emd_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND bessel_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND poisson_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND virial_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND average_bins_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND propermotion_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND EMD_Range_test_link_libs ${CRLIBM_LIBRARY})
    list(APPEND mixeddwarf_test_link_libs ${CRLIBM_LIBRARY})
endif()

milkyway_link(emd_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${emd_test_link_libs}")
milkyway_link(bessel_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${bessel_test_link_libs}")
milkyway_link(poisson_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${poisson_test_link_libs}")
milkyway_link(virial_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${virial_test_link_libs}")
milkyway_link(average_bins_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${average_bins_test_link_libs}")
milkyway_link(propermotion_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${propermotion_test_link_libs}")
milkyway_link(EMD_Range_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${EMD_Range_test_link_libs}")
milkyway_link(mixeddwarf_test ${BOINC_APPLICATION} ${NBODY_STATIC} "${mixeddwarf_test_link_libs}")

if(BOINC_APPLICATION)
  if(UNIX)
    target_link_libraries(nbody_test_driver pthread)
  endif()
endif()


add_test(NAME checkpoint_test
           WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests"
           COMMAND nbody_test_driver "CheckpointTest.lua")

add_test(NAME custom_arg_test
           WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests"
           COMMAND nbody_test_driver "RunArgumentTests.lua" $<TARGET_FILE:milkyway_nbody>)


add_test(NAME emd_test COMMAND emd_test)

add_test(NAME bessel_test COMMAND bessel_test)

add_test(NAME poisson_test COMMAND poisson_test)

add_test(NAME virial_test COMMAND virial_test)

add_test(NAME average_bins_test WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/test_hists" COMMAND average_bins_test)

add_test(NAME propermotion_test WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/test_hists" COMMAND propermotion_test)

add_test(NAME EMD_Range_test WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests/test_hists" COMMAND EMD_Range_test)

add_test(NAME mixeddwarf_test COMMAND mixeddwarf_test)

set(invalid_test_dir "${PROJECT_SOURCE_DIR}/tests/invalid_tests")
file(GLOB INVALID_TEST_INPUTS "${invalid_test_dir}/*.lua")
add_test(NAME invalid_input_test
           WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests"
           COMMAND nbody_test_driver "RunInvalidInputTests.lua"
                                   $<TARGET_FILE:milkyway_nbody>
                                   ${INVALID_TEST_INPUTS})

# Run these first 3 tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -I 1,4,1)

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
                       model_8
                       model_9
                       model_ninkovic
                       model_triaxial
                       model_newhist1
                       model_newhist2
                       model_newhist3
                       model_LMC
                       model_bar
                       model_LMC_bar)

foreach(model_name ${orphan_model_names})
  foreach(n ${body_counts})
    add_test(NAME ${model_name}__${n}_test
               WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests"
               COMMAND nbody_test_driver "RunTestUnits.lua"
                                          $<TARGET_FILE:milkyway_nbody>
                                          "orphan_models"
                                          ${model_name}
                                          "orphan_model_histogram_3"
                                          ${n})
  endforeach()
endforeach()

function(make_bodycount_test_set n)
  add_custom_target(test_${n} COMMAND ${CMAKE_CTEST_COMMAND} -R "model_.*__${n}_test$"
                              DEPENDS milkyway_nbody
                              COMMENT "Running model tests with ${n} bodies")
endfunction()

foreach(n ${body_counts})
  make_bodycount_test_set(${n})
endforeach()


#--------UNCOMMENT THIS SECTION TO ADD CONTEXT TEST--------#

add_test(NAME context_test
           WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests"
           COMMAND nbody_test_driver "ContextTest.lua")

#----------------------------------------------------------#

#compress: ${CMAKE_COMMAND} -E tar cvfj test_output.tar.bz2 test_output
#execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvfj "context_test_results.tar.bz2"
#                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")


add_custom_target(context_test COMMAND ${CMAKE_CTEST_COMMAND} -R "context_test")

add_custom_target(bench COMMAND nbody_test_driver "RunBenchmarks.lua"
                                                  $<TARGET_FILE:milkyway_nbody>
                                                  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests")


