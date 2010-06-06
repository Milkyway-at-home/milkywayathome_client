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
  set(BOINC_INCLUDE_SEARCH_PATH /usr/local/include/boinc )
endif(UNIX)

if(MINGW)
  set(BOINC_INCLUDE_SEARCH_PATH ${TAKEOFFGW_ROOT}/local/include
                                ${TAKEOFFGW_ROOT}/local/include/boinc)

  set(BOINC_LIB_SEARCH_PATH ${TAKEOFFGW_ROOT}/local/lib)
endif(MINGW)


find_path(BOINC_INCLUDE_DIR boinc_api.h ${BOINC_INCLUDE_SEARCH_PATH})
find_library(BOINC_LIBRARY boinc ${BOINC_LIB_SEARCH_PATH})
find_library(BOINC_API_LIBRARY boinc_api ${BOINC_LIB_SEARCH_PATH})

if(BOINC_LIBRARY AND BOINC_API_LIBRARY)
    set(BOINC_LIBRARIES ${BOINC_LIBRARY} ${BOINC_API_LIBRARY})
endif(BOINC_LIBRARY AND BOINC_API_LIBRARY)

if(BOINC_INCLUDE_DIR AND BOINC_LIBRARIES)
   set(BOINC_FOUND TRUE)
endif(BOINC_INCLUDE_DIR AND BOINC_LIBRARIES)


if(BOINC_FOUND)
   if(NOT Boinc_FIND_QUIETLY)
      message(STATUS "Found BOINC Libraries: ${BOINC_LIBRARIES}")
   endif(NOT Boinc_FIND_QUIETLY)
else(BOINC_FOUND)
   if(Boinc_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find BOINC Libraries")
   endif(Boinc_FIND_REQUIRED)
endif(BOINC_FOUND)


