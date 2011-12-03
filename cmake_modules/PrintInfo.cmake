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


macro(print_libs)
  message("   BOINC Libraries    ${BOINC_LIBRARIES}")
  message("   Lua Libraries      ${LUA_LIBRARIES}")
  message("   OpenCL Libraries   ${OPENCL_LIBRARIES}")
  message("   OpenGL Libraries   ${OPENGL_LIBRARIES}")
  message("   GLUT Libraries     ${GLUT_LIBRARY}")
  message("   OpenSSL Libraries  ${OPENSSL_LIBARIES}")
endmacro()

macro(print_build_options)
  message("   ----")
  message("   Double precision:  ${DOUBLEPREC}")
  message("   Separation:        ${SEPARATION}")
  message("   N-Body:            ${NBODY}")
  message("   OpenCL separation: ${SEPARATION_OPENCL}")
  message("   OpenCL nbody:      ${NBODY_OPENCL}")
  message("   BOINC application: ${BOINC_APPLICATION}")
  message("   BOINC release:     ${BOINC_RELEASE}")
  message("   ----")
endmacro()

macro(print_separator)
  message("---------------------------------------------")
endmacro()

macro(print_build_info)
  message("   Building:          ${PROJECT_NAME}")
  message("   System:            ${CMAKE_SYSTEM}")
  message("   Build type:        ${CMAKE_BUILD_TYPE}")
  message("   Arch:              ${CMAKE_SYSTEM_PROCESSOR}")
  if(APPLE)
    message("   OS X target:       ${CMAKE_OSX_DEPLOYMENT_TARGET}")
  endif()
  message("   Install path:      ${CMAKE_INSTALL_PREFIX}")
  message("   ----")
  message("   CMAKE version:     ${CMAKE_VERSION}")
  message("   CMAKE binary:      ${CMAKE_COMMAND}")
  message("   CTEST binary:      ${CMAKE_CTEST_COMMAND}")
  message("   CMAKE generator:   ${CMAKE_GENERATOR}")
  message("   ----")
  message("   Project src dir:   ${CMAKE_SOURCE_DIR}")
  message("   Project bin dir:   ${CMAKE_BINARY_DIR}")
  message("   Build tool:        ${CMAKE_BUILD_TOOL}")
  message("   C Compiler:        ${CMAKE_C_COMPILER}")
  #TODO:Report CFLAGS used based on build type
  #message("   CFLAGS:           ${CMAKE_C_FLAGS}")
endmacro()

