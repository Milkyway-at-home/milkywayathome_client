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

if(CPU_SEPARATION)
  add_executable(milkyway_client milkyway_client_app.c)
  set_target_properties(milkyway_client
                          PROPERTIES
                            COMPILE_FLAGS "-DCOMPUTE_ON_CPU=1 -DUSE_CUDA=0 -D USE_OCL=0")
  target_link_libraries(milkyway_client ${POPT_LIBRARY} ${BOINC_LIBRARIES} milkyway_separation_cpu)
  install(TARGETS milkyway_client DESTINATION bin
       PERMISSIONS
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)
endif()


if(OPENCL_SEPARATION)
  add_executable(milkyway_client_ocl milkyway_client_app.c)
  target_link_libraries(milkyway_client_ocl ${POPT_LIBRARY}
                                            ${BOINC_LIBRARIES}
                                            milkyway_separation_ocl
                                            milkyway_separation_cpu)
  set_target_properties(milkyway_client_ocl
                        PROPERTIES
                        COMPILE_FLAGS "-DCOMPUTE_ON_CPU=0 -DUSE_CUDA=0 -D USE_OCL=1")

  install(TARGETS milkyway_client_ocl DESTINATION bin
       PERMISSIONS
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)
endif()


if(CUDA_SEPARATION)
  add_executable(milkyway_client_cuda milkyway_client_app.c)
  target_link_libraries(milkyway_client_ocuda${POPT_LIBRARY}
                                            ${BOINC_LIBRARIES}
                                            milkyway_separation_cuda
                                            milkyway_separation_cpu)
  set_target_properties(milkyway_client_cuda
                        PROPERTIES
                        COMPILE_FLAGS "-DCOMPUTE_ON_CPU=0 -DUSE_CUDA=1 -D USE_OCL=0")
  install(TARGETS milkyway_client_cuda DESTINATION bin
       PERMISSIONS
       OWNER_READ OWNER_WRITE OWNER_EXECUTE
       GROUP_READ GROUP_EXECUTE
       WORLD_READ WORLD_EXECUTE)
endif()


