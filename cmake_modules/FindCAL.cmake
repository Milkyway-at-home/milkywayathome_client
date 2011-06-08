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

find_library(ATICALCL_LIBRARY aticalcl PATHS "$ENV{AMDAPPSDKROOT}/lib")
find_library(ATICALRT_LIBRARY aticalrt PATHS "$ENV{AMDAPPSDKROOT}/lib")
find_path(ATICAL_INCLUDE_DIR NAMES CAL/cal.h CAL/calcl.h PATHS "$ENV{AMDAPPSDKROOT}/include")

if(ATICAL_INCLUDE_DIR AND ATICALCL_LIBRARY AND ATICALRT_LIBRARY)
   set(CAL_FOUND TRUE)
endif()

if(CAL_FOUND)
   if(NOT CAL_FIND_QUIETLY)
      message(STATUS "Found calcl Library: ${ATICALCL_LIBRARY}")
      message(STATUS "Found calrt Library: ${ATICALRT_LIBRARY}")
   endif()
else()
   if(CAL_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find CAL Library")
   endif()
endif()

