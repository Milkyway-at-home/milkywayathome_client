# Copyright 2011 Matthew Arsenault, Travis Desell, Dave Przybylo,
# Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
# Magdon-Ismail and Rensselaer Polytechnic Institute.
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


if("${CMAKE_C_COMPILER_ID}" MATCHES "PathScale")
  set(C_COMPILER_IS_PATHCC 1)
endif()

if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  set(C_COMPILER_IS_CLANG 1)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CXX_COMPILER_IS_CLANG 1)
endif()

if(CXX_COMPILER_IS_CLANG OR C_COMPILER_IS_CLANG)
  set(COMPILER_IS_CLANG 1)
endif()


if("${CMAKE_C_COMPILER_ID}" MATCHES "SunPro")
  set(C_COMPILER_IS_SUN 1)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "SunPro")
  set(CXX_COMPILER_IS_SUN 1)
endif()


