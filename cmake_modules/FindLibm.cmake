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


if(NOT MSVC)
  find_path(LIBM_INCLUDE_DIR math.h)
else()
  set(LIBM_INCLUDE_DIR ${CMAKE_INCLUDE_PATH})
endif()

if(APPLE)
  find_library(LIBM_LIBRARY m)
elseif(WIN32)
  set(LIBM_LIBRARY "")
else()
  if(LIBM_USE_STATIC)
    set(__old_cmake_find_lib_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()

  find_library(LIBM_LIBRARY m)
  if(LIBM_USE_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${__old_cmake_find_lib_suffixes})
  endif()

  if(LIBM_INCLUDE_DIR AND LIBM_LIBRARY)
    set(LIBM_FOUND TRUE)
  endif()

  if(LIBM_FOUND)
    if(NOT Libm_FIND_QUIETLY)
      message(STATUS "Found LIBM Library: ${LIBM_LIBRARY}")
    endif(NOT Libm_FIND_QUIETLY)
  else()
    if(Libm_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find LIBM Library")
    endif(Libm_FIND_REQUIRED)
  endif()
endif()

