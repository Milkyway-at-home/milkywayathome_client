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

option(BUILD_SHARED_LIBS
       "Set to OFF to build static libraries" ON)

option(DOUBLE_PRECISION
       "Use double precision" ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(UNIX AND NOT APPLE)  #Unix that isn't OS X
    set(OS_SPECIFIC_LIBS pthread)
endif(UNIX AND NOT APPLE)

if(WIN32)
  set(OS_SPECIFIC_LIBS msvcrt)
endif(WIN32)

set(lib_source_file_list evaluation_state.c
                         simple_evaluator.c
                         ${CMAKE_SOURCE_DIR}/util/io_util.c
                         atSurveyGeometry.c
                         stCoords.c
                         stVector.c
                         numericalIntegration.c
                         evaluation_optimized.c
                         parameters.c
                         star_points.c
                         astronomy_worker.c)

# Windows is a pain in the ass.
if(WIN32)
  add_library(milkyway_separation_cpu STATIC ${lib_source_file_list})
else(WIN32)
  add_library(milkyway_separation_cpu SHARED ${lib_source_file_list})
endif(WIN32)

target_link_libraries(milkyway_separation_cpu ${BOINC_LIBRARIES} ${OS_SPECIFIC_LIBS} m)

install(TARGETS milkyway_separation_cpu DESTINATION lib
       PERMISSIONS           #Fix occasional issuse with people who change their umask
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_INCLUDE_DIR}/milkyway.h DESTINATION include/milkyway)
install(FILES ${PROJECT_INCLUDE_DIR}/config.h DESTINATION include/milkyway)

