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

include(CheckCCompilerFlag)

function(add_flag_if_supported flagname)
  message("Checking for flag ${flagname}")

  check_c_compiler_flag("${flagname}" HAVE_FLAG_${flagname})

  if(${HAVE_FLAG_${flagname}})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flagname}" PARENT_SCOPE)
  endif()
endfunction()

function(append_supported_flags flags)
  foreach(flag ${flags})
    add_flag_if_supported(${flag})
  endforeach()
endfunction()

