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

cmake_minimum_required (VERSION 2.6)

set(test_app ${PROJECT_BINARY_DIR}/bin/milkyway_client)

#arguments shared by several of the tests
set(param_set1 0.571713 12.312119 -3.305187 148.010257 22.453902 0.42035
               -0.468858 0.760579 -1.361644 177.884238 23.882892 1.210639
               -1.611974 8.534378 -1.361644 177.884238 10.882892 1.210639
               -1.611974 8.534378)

set(test_11_params ${param_set1})
set(test_12_params ${param_set1})
set(test_20_params ${param_set1})
set(test_21_params ${param_set1})

set(test_79_params 0.34217373320392042 25.9517910846623 -2.1709414738826602 38.272511356953906
                   30.225190442596112 2.2149060013372885 0.32316169064291655 2.7740244716285285)

set(test_82_params 0.40587961154742185 17.529961843393409 -1.8575145272144837 29.360893891378243
                   31.228263575178566 -1.551741065334 .064096152599308373 2.55428209912781)

set(test_86_params 0.73317163557524425 14.657212876628332 -1.7054653473950408 16.911711745343633
                   28.077212666463502 -1.2032908515814611 3.5273606439247287 2.2248214505875008)

set(test_list 11 12 20 21 79 82 86)

foreach(i ${test_list})
  if(SMALL_TESTS)
    set(small_test_${i}_args -p stars-${i}.txt -a astronomy_parameters-${i}-small.txt ${test_${i}_params})
    add_test(small_test_${i} ${test_app} ${small_test_${i}_args})
  endif()

  if(LARGE_TESTS)
    set(test_${i}_args -p stars-${i}.txt -a astronomy_parameters-${i}.txt ${test_${i}_params})
    add_test(test_${i} ${test_app} ${test_${i}_args})
  endif()

  if(AUX_SMALL_TESTS)
    set(small_test_${i}_args-aux -p stars-${i}.txt -a astronomy_parameters-${i}-small-aux.txt ${test_${i}_params})
    add_test(small_test_${i}-aux ${test_app} ${small_test_${i}_args-aux})
  endif()

  if(AUX_LARGE_TESTS)
    set(test_${i}_args-aux -p stars-${i}.txt -a astronomy_parameters-${i}-aux.txt ${test_${i}_params})
    add_test(test_${i}-aux ${test_app} ${test_${i}_args-aux})
  endif()

endforeach()



# List of tests to run.
# Test sources should be named <testname>.c .
set(test_list )


foreach(test ${test_list})
	add_executable(${test} ${test}.c ${extra_test_sources})
	add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
    add_dependencies(check ${test})
    target_link_libraries(${test} m)
endforeach()

