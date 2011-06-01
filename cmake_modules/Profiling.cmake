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


# set(CMAKE_BUILD_TYPE Profile)

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS_PROFILE "-g -pg" CACHE STRING "C compiler flags for profiling builds" FORCE)
  set(CMAKE_CXX_FLAGS_PROFILE "-g -pg" CACHE STRING "C++ compiler flags for profiling builds" FORCE)

  set(CMAKE_EXE_LINKER_FLAGS_PROFILE "" CACHE STRING "" FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "" CACHE STRING "" FORCE)
else()
  # TODO
  set(CMAKE_C_FLAGS_PROFILE CACHE STRING "C compiler flags for profiling builds")
  set(CMAKE_CXX_FLAGS_PROFILE CACHE STRING "C++ compiler flags for profiling builds")

  set(CMAKE_EXE_LINKER_FLAGS_PROFILE "" CACHE STRING "" FORCE)
  set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "" CACHE STRING "" FORCE)
endif()

mark_as_advanced(CMAKE_C_FLAGS_PROFILE
                 CMAKE_CXX_FLAGS_PROFILE
                 CMAKE_EXE_LINKER_FLAGS_PROFILE
                 CMAKE_SHARED_LINKER_FLAGS_PROFILE)

