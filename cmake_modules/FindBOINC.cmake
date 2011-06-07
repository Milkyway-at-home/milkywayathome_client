# Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
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

include(CheckStructHasMember)

if(UNIX)
  set(BOINC_INCLUDE_SEARCH_PATH /usr/local/include/boinc
                                /usr/local/include
                                /usr/include/boinc)
endif(UNIX)

if(MINGW)
  set(BOINC_INCLUDE_SEARCH_PATH /local/include
                                /local/include/boinc)

  set(BOINC_LIB_SEARCH_PATH /local/lib)
endif(MINGW)

find_path(BOINC_INCLUDE_DIR boinc_api.h ${BOINC_INCLUDE_SEARCH_PATH})

if(NOT MSVC)
  find_library(BOINC_LIBRARY boinc ${BOINC_LIB_SEARCH_PATH})
  find_library(BOINC_API_LIBRARY boinc_api ${BOINC_LIB_SEARCH_PATH})
  find_library(BOINC_GRAPHICS_LIBRARY boinc_graphics2 ${BOINC_LIB_SEARCH_PATH})
else()
  if(NOT MSVC_USE_STATIC_CRT)
    find_library(BOINC_LIBRARY libboinc ${BOINC_LIB_SEARCH_PATH})
    find_library(BOINC_API_LIBRARY libboincapi ${BOINC_LIB_SEARCH_PATH})
    find_library(BOINC_GRAPHICS_LIBRARY boinc_graphics2 ${BOINC_LIB_SEARCH_PATH})
  else()
    find_library(BOINC_LIBRARY libboinc_staticcrt ${BOINC_LIB_SEARCH_PATH})
    find_library(BOINC_API_LIBRARY libboincapi_staticcrt ${BOINC_LIB_SEARCH_PATH})
    find_library(BOINC_GRAPHICS_LIBRARY boinc_graphics2 ${BOINC_LIB_SEARCH_PATH})
  endif()
endif()

set(BOINC_LIBRARIES ${BOINC_LIBRARY})

if(BOINC_API_LIBRARY)
  list(INSERT BOINC_LIBRARIES 0 ${BOINC_API_LIBRARY})
endif()

if(BOINC_GRAPHICS_LIBRARY)
  list(INSERT BOINC_LIBRARIES 0 ${BOINC_GRAPHICS_LIBRARY})
endif()

if(BOINC_GRAPHICS_LIBRARY)
  set(BOINC_GRAPHICS_FOUND TRUE)
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

set(CMAKE_REQUIRED_INCLUDES ${BOINC_INCLUDE_DIR})
check_struct_has_member("BOINC_OPTIONS" "multi_thread" "boinc_api.h" HAVE_BOINC_OPTIONS_MULTI_THREAD)
if(NOT HAVE_BOINC_OPTIONS_MULTI_THREAD)
  message(FATAL_ERROR "Found BOINC libraries too old")
endif()
set(CMAKE_REQUIRED_INCLUDES)

