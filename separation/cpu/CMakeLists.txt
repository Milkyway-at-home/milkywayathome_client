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

set(lib_source_file_list evaluation_state.c
                         ${PROJECT_SOURCE_DIR}/util/io_util.c
                         coordinates.c
                         numerical_integration.c
                         evaluation.c
                         parameters.c
                         star_points.c
                         integrals.c
                         likelihood.c
                         integral_constants.c
                         separation_cl_buffers.c
                         setup_cl.c
                         run_cl.c
                         separation_cl_defs.c)

add_library(separation_cpu STATIC ${lib_source_file_list})

include_directories(${OPENCL_INCLUDE_DIRS})


target_link_libraries(separation_cpu ${BOINC_LIBRARIES} ${OS_SPECIFIC_LIBS} m ${OPENCL_LIBRARIES})

install(TARGETS separation_cpu DESTINATION lib
       PERMISSIONS           #Fix occasional issuse with people who change their umask
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)

install(FILES ${PROJECT_INCLUDE_DIR}/separation.h DESTINATION include/milkyway)
install(FILES ${PROJECT_INCLUDE_DIR}/config.h DESTINATION include/milkyway)

