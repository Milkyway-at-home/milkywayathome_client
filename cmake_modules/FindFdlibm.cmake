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


find_path(FDLIBM_INCLUDE_DIR fdlibm.h)
find_library(FDLIBM_LIBRARY fdlibm)

if(FDLIBM_INCLUDE_DIR AND FDLIBM_LIBRARY)
  set(FDLIBM_FOUND TRUE)
endif()

if(FDLIBM_FOUND)
  if(NOT Fdlibm_FIND_QUIETLY)
    message(STATUS "Found Fdlibm Library: ${FDLIBM_LIBRARY}")
  endif(NOT Fdlibm_FIND_QUIETLY)
else()
  if(Fdlibm_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find Fdlibm Library")
  endif(Fdlibm_FIND_REQUIRED)
endif()


