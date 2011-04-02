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
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.6)
project(lmodl)

enable_testing()


find_package(SDL REQUIRED)

include_directories(include ${SDL_INCLUDE_DIR})

set(LMODL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" CACHE STRING "LModL includes")
include_directories(${LMODL_INCLUDE_DIR})


message("---------------------------------------------")
message("   Building:          ${PROJECT_NAME}")
message("   System:            ${CMAKE_SYSTEM}")
message("   Build type:        ${CMAKE_BUILD_TYPE}")
message("   Arch:              ${CMAKE_SYSTEM_PROCESSOR}")
message("   Install path:      ${CMAKE_INSTALL_PREFIX}")
message("   ----")
message("   CMAKE version:     ${CMAKE_VERSION}")
message("   CMAKE binary:      ${CMAKE_COMMAND}")
message("   CTEST binary:      ${CMAKE_CTEST_COMMAND}")
message("   CMAKE generator:   ${CMAKE_GENERATOR}")
message("   ----")
message("   Project src dir:   ${CMAKE_SOURCE_DIR}")
message("   Project bin dir:   ${CMAKE_BINARY_DIR}")
message("   ----")
message("   SDL Library        ${SDL_LIBRARY}")
message("   ----")
message("   Build tool:        ${CMAKE_BUILD_TOOL}")
message("   C Compiler:        ${CMAKE_C_COMPILER}")
message("   C++ Compiler       ${CMAKE_CXX_COMPILER}")
message("---------------------------------------------")

set(lmodl_headers include/3dproj.hpp
                  include/draw.hpp
                  include/draw32.hpp
                  include/draw8.hpp
                  include/drawhalo.hpp
                  include/quatmath.hpp
                  include/trigtbl.hpp
                  include/binfile.hpp
                  include/draw24.hpp
                  include/draw32p.hpp
                  include/drawcore.hpp
                  include/imgplot.hpp
                  include/rgbconv.hpp
                  include/vectmath.hpp)

set(lmodl_src src/3dproj.cpp
              src/draw.cpp
              src/draw32p.cpp
              src/drawcore.cpp
              src/imgplot.cpp
              src/trigtbl.cpp
              src/binfile.cpp
              src/draw32.cpp
              src/draw8.cpp
              src/drawhalo.cpp
              src/rgbconv.cpp)

add_library(lmodl STATIC ${lmodl_src} ${lmodl_headers})

set(test_list fpstest
              gfxinfo
              imgrender
              quattest
              trigtblm)


foreach(i ${test_list})
  add_executable(${i} test/${i}.cpp)
  target_link_libraries(${i} ${SDL_LIBRARY} lmodl)
endforeach()


