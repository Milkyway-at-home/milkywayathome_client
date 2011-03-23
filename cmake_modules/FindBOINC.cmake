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

#FIXME: This isn't really good.
if(UNIX)
  set(BOINC_INCLUDE_SEARCH_PATH /usr/local/include/boinc
                                /usr/include/boinc)
endif(UNIX)

if(MINGW)
  set(BOINC_INCLUDE_SEARCH_PATH ${TAKEOFFGW_ROOT}/local/include
                                ${TAKEOFFGW_ROOT}/local/include/boinc)

  set(BOINC_LIB_SEARCH_PATH ${TAKEOFFGW_ROOT}/local/lib)
endif(MINGW)

if(BOINC_USE_STATIC)
  set(__old_cmake_find_lib_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

find_path(BOINC_INCLUDE_DIR boinc/boinc_api.h ${BOINC_INCLUDE_SEARCH_PATH})

if(NOT MSVC)
  find_library(BOINC_LIBRARY boinc ${BOINC_LIB_SEARCH_PATH})
  find_library(BOINC_API_LIBRARY boinc_api ${BOINC_LIB_SEARCH_PATH})
else()
  if(NOT MSVC_USE_STATIC_CRT)
    find_library(BOINC_LIBRARY libboinc ${BOINC_LIB_SEARCH_PATH})
    find_library(BOINC_API_LIBRARY libboincapi ${BOINC_LIB_SEARCH_PATH})
  else()
    find_library(BOINC_LIBRARY libboinc_staticcrt ${BOINC_LIB_SEARCH_PATH})
    find_library(BOINC_API_LIBRARY libboincapi_staticcrt ${BOINC_LIB_SEARCH_PATH})
  endif()
endif()

if(BOINC_USE_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${__old_cmake_find_lib_suffixes})
endif()

set(BOINC_LIBRARIES ${BOINC_LIBRARY})

if(BOINC_API_LIBRARY)
  list(INSERT BOINC_LIBRARIES 0 ${BOINC_API_LIBRARY})
endif()


if(BOINC_INCLUDE_DIR AND BOINC_LIBRARY)
   set(BOINC_FOUND TRUE)
endif()


if(BOINC_FOUND)
   if(NOT Boinc_FIND_QUIETLY)
      message(STATUS "Found BOINC Libraries: ${BOINC_LIBRARY}")
   endif(NOT Boinc_FIND_QUIETLY)
else(BOINC_FOUND)
   if(Boinc_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find BOINC Libraries")
   endif(Boinc_FIND_REQUIRED)
endif()


