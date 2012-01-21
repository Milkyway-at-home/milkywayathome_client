#
# Copyright 2011 Matthew Arsenault, Rensselaer Polytechnic Institute.
#
# This file is part of Milkway@Home.
#
# Milkyway@Home is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Milkyway@Home is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.8)
project(separation_il_kernelgen)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")



find_package(Boost)
find_package(CAL)
find_package(CAL++)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${ATICAL_INCLUDE_DIR})
include_directories(${ATICAL_INCLUDE_DIR}/CAL)
include_directories(${CAL++_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
add_executable(il_kernelgen il_kernelgen.cpp)
target_link_libraries(il_kernelgen ${ATICALRT_LIBRARY})


set(MILKYWAYATHOME_CLIENT_CMAKE_MODULES "${milkywayathome_client_SOURCE_DIR}/cmake_modules")
list(APPEND CMAKE_MODULE_PATH ${MILKYWAYATHOME_CLIENT_CMAKE_MODULES})


