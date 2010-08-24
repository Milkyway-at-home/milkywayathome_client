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

cmake_minimum_required(VERSION 2.8)
project(milkyway_separation)

enable_testing()

set(milkyway_separation_VERSION_MAJOR 0)
set(milkyway_separation_VERSION_MINOR 4)
set(milkyway_separation_${milkyway_client_MAJOR_VERSION}.${milkyway_client_MINOR_VERSION})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include_directories(${PROJECT_INCLUDE_DIR})
set(CMAKE_C_FLAGS "-std=c99 -warn-uninit-values ${CMAKE_C_FLAGS}")

configure_file (
  "${PROJECT_INCLUDE_DIR}/config.h.in"
  "${PROJECT_INCLUDE_DIR}/config.h")


if(UNIX AND NOT APPLE)  #Unix that isn't OS X
    set(OS_SPECIFIC_LIBS pthread)
endif()

if(WIN32)
  set(OS_SPECIFIC_LIBS msvcrt)
endif()

add_definitions("-DSEPARATION_OPENCL=0")
add_definitions("-DDOUBLEPREC=1")


add_subdirectory(tests EXCLUDE_FROM_ALL)


set(lib_source_file_list src/evaluation_state.c
                         src/io_util.c
                         src/gauss_legendre.c
                         src/evaluation.c
                         src/parameters.c
                         src/star_points.c
                         src/integrals.c
                         src/likelihood.c
                         src/calculated_constants.c
                         src/r_points.c
                         src/separation_cl_buffers.c
                         src/setup_cl.c
                         src/run_cl.c
                         src/separation_cl_defs.c)

add_library(separation_cpu STATIC ${lib_source_file_list})

include_directories(${OPENCL_INCLUDE_DIRS})


#target_link_libraries(separation_cpu ${BOINC_LIBRARIES} ${OS_SPECIFIC_LIBS} m ${OPENCL_LIBRARIES})


# install(TARGETS separation_cpu DESTINATION lib
#        PERMISSIONS           #Fix occasional issuse with people who change their umask
#        OWNER_READ OWNER_WRITE OWNER_EXECUTE
#        GROUP_READ GROUP_EXECUTE
#        WORLD_READ WORLD_EXECUTE)

# install(FILES ${PROJECT_INCLUDE_DIR}/separation.h DESTINATION include/milkyway)
# install(FILES ${PROJECT_INCLUDE_DIR}/config.h DESTINATION include/milkyway)


set(client_bin_name milkyway_separation)

if(NEED_POPT_DEPS)
  # Popt dependencies for static linking For some reason, on OS X /
  # windows you need these dependencies but not Linux.
  list(APPEND exe_link_libs ${LIBICONV_LIBRARY}
                            ${LIBINTL_LIBRARY})
endif()



add_executable(${client_bin_name} src/milkyway_client_app.c)
set_target_properties(${client_bin_name}
                      PROPERTIES
                       COMPILE_FLAGS "-DCOMPUTE_ON_CPU=1 -DUSE_CUDA=0 -D USE_OCL=0")
  target_link_libraries(${client_bin_name} ${POPT_LIBRARY} ${BOINC_LIBRARIES} separation_cpu milkyway
                           ${exe_link_libs})

  install(TARGETS ${client_bin_name} DESTINATION bin
          PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE)


  if(NOT APPLE)
    #CHECKME: What about Windows?
    set(separation_static_link_flags "-static -static-libgcc -static-libstdc++")

    if(UNIX)
      #We need to do this to statically link pthreads. Otherwise really
      #bizarre things happen.
      set(separation_static_link_flags "-pthread ${nbody_static_link_flags}")
    endif()

    set_target_properties(${client_bin_name}
      PROPERTIES
      LINKER_LANGUAGE CXX
      LINK_SEARCH_END_STATIC ON
      LINK_FLAGS ${separation_static_link_flags})
  else()
    set_target_properties(${client_bin_name}
      PROPERTIES
      LINKER_LANGUAGE CXX
      LINK_SEARCH_END_STATIC ON)
  endif()




if(OPENCL_SEPARATION)
  add_executable(${client_bin_name}_ocl milkyway_client_app.c)
  target_link_libraries(${client_bin_name}_ocl ${POPT_LIBRARY}
                                               ${BOINC_LIBRARIES}
                                               milkyway_separation_ocl
                                               milkyway_separation_cpu)
  set_target_properties(${client_bin_name}_ocl
                        PROPERTIES
                        COMPILE_FLAGS "-DCOMPUTE_ON_CPU=0 -DUSE_CUDA=0 -D USE_OCL=1")

  install(TARGETS ${client_bin_name}_ocl DESTINATION bin
          PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE
            WORLD_READ WORLD_EXECUTE)
endif()


