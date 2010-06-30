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

cmake_minimum_required (VERSION 2.8)
project(milkyway_nbody)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


if(CMAKE_COMPILER_IS_GNUCC)
  set(DEFAULT_USE_DYNAMIC_PRECISION ON)
else()
  set(DEFAULT_USE_DYNAMIC_PRECISION OFF)
endif()

option(DYNAMIC_PRECISION
       "Use dynamic precision" ${DEFAULT_USE_DYNAMIC_PRECISION})

option(BOINC_APPLICATION
  "Build as a BOINC Application" ON)

option(BOINC_DEBUG
  "Avoid annoying things BOINC does for debugging" ON)

#Used if not using dynamic precision
option(DOUBLEPREC
       "Use double precision" ON)

option(OUTPUT_CARTESIAN
       "Output galactic Cartesian coordinates instead of lbr" OFF)

if(OUTPUT_CARTESIAN)
  add_definitions("-DOUTPUT_CARTESIAN")
endif()

set(BOINC_DEFS "")

if(BOINC_APPLICATION)
  list(APPEND BOINC_DEFS "BOINC_APPLICATION=1")
else()
  list(APPEND BOINC_DEFS "BOINC_APPLICATION=0")
endif()

if(BOINC_DEBUG)
  list(APPEND BOINC_DEFS "BOINC_DEBUG=1")
else()
  list(APPEND BOINC_DEFS "BOINC_DEBUG=0")
endif()


set(nbody_lib_src ${PROJECT_SOURCE_DIR}/src/chisq.c
                  ${PROJECT_SOURCE_DIR}/src/grav.c
                  ${PROJECT_SOURCE_DIR}/src/io.c
                  ${PROJECT_SOURCE_DIR}/src/load.c
                  ${PROJECT_SOURCE_DIR}/src/orbitintegrator.c
                  ${PROJECT_SOURCE_DIR}/src/nbody_simulation.c
                  ${PROJECT_SOURCE_DIR}/src/json_params.c
                  ${PROJECT_SOURCE_DIR}/src/plummer.c
                  ${PROJECT_SOURCE_DIR}/src/nbody_util.c)

set(CMAKE_C_FLAGS "-std=gnu99 -fvisibility=hidden ${CMAKE_C_FLAGS}")



message("---------------------------------------------")
message("   Building:          ${PROJECT_NAME}")
message("   Dynamic precision: ${DYNAMIC_PRECISION}")
message("   Double precision:  ${DOUBLEPREC}")
message("---------------------------------------------")

if(UNIX AND NOT APPLE)  #Unix that isn't OS X
    set(OS_SPECIFIC_LIBS pthread)
endif()

set(nbody_link_libs ${JSON_C_LIBRARY} m ${OS_SPECIFIC_LIBS})

if(BOINC_APPLICATION)
  list(APPEND nbody_link_libs ${BOINC_LIBRARIES})
endif()

# This dynamic precision is kind of sketchy. We compile the library
# twice, once with float and once with double precision. Anything
# using real shouldn't be exposed to the main application. The only
# symbol which (ideally) would be visible from the libraries would be
# the runNBodySimulation_* functions. All the other function names are
# the same, since this is actually the same library. This mostly works
# fine with a shared library, since symbols are resolved within the
# library when you build it, and functions within the library will
# link to the other functions internal to the right version of the
# library when it's compiled.

# It turns out static libraries are dumber than I thought, and don't
# do any sort of linking within themselves when you build them for
# lack of a better term. This means that when you try to link with the
# 2 static libraries, the first one you use resolves all the symbols
# in the first library. For the next one, since all the same symobls
# have been found in the first version, the only symbol left is the
# specific type runNBodySimulation_whatever. Then the second library
# will try to use the functions from the first which we would like to
# be private, and not the ones expecting the same type within the
# library. The end result is a different type of failure depending on
# whether you try linking the float or double version first. If you
# try using the double version first, the double will work. When you
# try to use single precision, due to the wrong functions being used
# with the smaller size, you end up with corruption in everything. If
# you specify them in the other order (the float library first, then
# the double), you get segfaults from trying to use structs that are
# bigger than the functions the linker found are expecting.

#TLDR: to use static libraries, need to figure out how to use some
#kind of symbol visibility trickery to expose only the necessary
#runNBodySimulation* symbols, as a shared library does for you.
#something involving -fvisibility=hidden and different sets of .o
#files. Dynamic precision will likely only be useful for testing
#things, and we won't bother with it for releases. It approximately
#doubles the binary size anyway.
if(DYNAMIC_PRECISION)
  add_definitions("-DDYNAMIC_PRECISION=1")

  add_library(nbody_simulation_double SHARED ${nbody_lib_src})
  target_link_libraries(nbody_simulation_double ${nbody_link_libs})
  set_target_properties(nbody_simulation_double
                          PROPERTIES
                            COMPILE_DEFINITIONS "DOUBLEPREC;${BOINC_DEFS}")


  add_library(nbody_simulation_float SHARED ${nbody_lib_src})
  target_link_libraries(nbody_simulation_float ${nbody_link_libs})

  set_target_properties(nbody_simulation_float
                          PROPERTIES
                            COMPILE_DEFINITIONS "${BOINC_DEFS}")


  set(nbody_libs nbody_simulation_float nbody_simulation_double)
else()
  if(DOUBLEPREC)
    add_definitions("-DDOUBLEPREC")
  endif()

  add_library(nbody_simulation STATIC ${nbody_lib_src})

  set_target_properties(nbody_simulation
                          PROPERTIES
                          COMPILE_DEFINITIONS "${BOINC_DEFS}")


  set(nbody_libs nbody_simulation)
endif()

list(APPEND nbody_link_libs ${nbody_libs})



add_executable(milkyway_nbody ${PROJECT_SOURCE_DIR}/src/main.c)
target_link_libraries(milkyway_nbody ${nbody_link_libs} ${POPT_LIBRARY})

set_target_properties(milkyway_nbody
                        PROPERTIES
                        COMPILE_DEFINITIONS "${BOINC_DEFS}")

add_subdirectory(tests EXCLUDE_FROM_ALL)

