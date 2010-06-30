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

foreach(test ${test_list})
	add_executable(${test} ${test}.c ${nbody_lib_src})
	add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
    add_dependencies(check ${test})
    target_link_libraries(${test} ${nbody_link_libs})
endforeach()


#Run the whole simulation with the test inputs, and compare the results.
#Use both double and single precision if available.

set(nbody_bin ${EXECUTABLE_OUTPUT_PATH}/milkyway_nbody)

set(result_dir /tmp/test_results)
set(nbody_test_dir ${CMAKE_CURRENT_SOURCE_DIR})

function(nbody_sample_test test_set test_set_name)
  set(double_out_dir ${result_dir}/${test_set_name}/double)
  set(float_out_dir ${result_dir}/${test_set_name}/float)
  file(MAKE_DIRECTORY ${double_out_dir} ${float_out_dir})

  foreach(i ${test_set})
    get_filename_component(fname ${i} NAME)
    get_filename_component(testname ${fname} NAME_WE)

    function(add_double_tests arg)
      set(outfile_double ${double_out_dir}/${testname}_results)
      set(reference_double ${nbody_test_dir}/${test_set_name}/double/${testname}_results)

      add_test(double_${testname} ${nbody_bin} ${arg} -f ${i} -o ${outfile_double})
      add_test(compare_double_${testname} ${CMAKE_COMMAND} -E compare_files ${reference_double} ${outfile_double})
    endfunction()

    function(add_float_tests)
      set(outfile_float ${float_out_dir}/${testname}_results)
      set(reference_float  ${nbody_test_dir}/${test_set_name}/float/${testname}_results)

      add_test(float_${testname} ${nbody_bin} -f ${i} -o ${outfile_float})
      add_test(compare_float_${testname} ${CMAKE_COMMAND} -E compare_files ${reference_float} ${outfile_float})
    endfunction()

    if(DYNAMIC_PRECISION)
      add_double_tests("-d")
      add_float_tests()
    else()
      if(DOUBLEPREC)
        add_double_tests("")
      else()
        add_float_tests()
      endif()
    endif()

  endforeach()
endfunction()

file(GLOB test_set1 "test_set1/*.js")
file(GLOB test_set2 "test_set2/*.js")

nbody_sample_test("${test_set1}" "test_set1")
nbody_sample_test("${test_set2}" "test_set2")

