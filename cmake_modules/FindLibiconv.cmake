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


find_path(LIBINTL_INCLUDE_DIR iconv.h)

find_library(LIBICONV_LIBRARY iconv)

if(LIBICONV_INCLUDE_DIR AND LIBICONV_LIBRARY)
   set(LIBICONV_FOUND TRUE)
endif(LIBICONV_INCLUDE_DIR AND LIBICONV_LIBRARY)

if(LIBICONV_FOUND)
   if(NOT Iconv_FIND_QUIETLY)
      message(STATUS "Found LIBICONV Library: ${LIBICONV_LIBRARY}")
   endif(NOT Iconv_FIND_QUIETLY)
else(LIBICONV_FOUND)
   if(Iconv_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find LIBICONV Library")
   endif(Iconv_FIND_REQUIRED)
endif(LIBICONV_FOUND)

