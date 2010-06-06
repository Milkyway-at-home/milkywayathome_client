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


find_path(POPT_INCLUDE_DIR popt.h)
find_library(POPT_LIBRARY popt)

if(POPT_INCLUDE_DIR AND POPT_LIBRARY)
   set(POPT_FOUND TRUE)
endif(POPT_INCLUDE_DIR AND POPT_LIBRARY)

if(POPT_FOUND)
   if(NOT Popt_FIND_QUIETLY)
      message(STATUS "Found POPT Library: ${POPT_LIBRARY}")
   endif(NOT Popt_FIND_QUIETLY)
else(POPT_FOUND)
   if(Popt_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find POPT Library")
   endif(Popt_FIND_REQUIRED)
endif(POPT_FOUND)

