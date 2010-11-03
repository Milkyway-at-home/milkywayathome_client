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

cmake_minimum_required(VERSION 2.8)
project(simulation)

find_package(SDL REQUIRED)
find_package(SDL_image REQUIRED)

include_directories(${SDL_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/LModL/include)
link_directories(${PROJECT_SOURCE_DIR}/LModL/lib)


add_subdirectory(LModL)





set(test_list mwdemo_fullscreen)


add_executable(nbody_demo bin/nbody.cpp bin/mwdemo.cpp)
target_link_libraries(nbody_demo lmodl)

# foreach(i ${test_list})
#   add_executable(${i} bin/${i}.cpp)
#   target_link_libraries(${i} lmodl)
# endforeach()





