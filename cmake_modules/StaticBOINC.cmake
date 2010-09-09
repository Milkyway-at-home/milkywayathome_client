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

# For releases, we want to be able to statically link as much as
# possible. This requires special handling on OS X since Apple
# doesn't like you trying to statically link the standard libraries.
# We also have to link as C++ when we do this because of BOINC.

function(correct_static_link client_bin_name)
  if(NOT APPLE)
    #CHECKME: What about Windows?
    set(client_static_link_flags "-static -static-libgcc -static-libstdc++")

    if(UNIX)
      #We need to do this to statically link pthreads. Otherwise really
      #bizarre things happen.
      set(client_static_link_flags "-pthread ${client_static_link_flags}")
    endif()

    set_target_properties(${client_bin_name}
      PROPERTIES
      LINKER_LANGUAGE CXX
      LINK_SEARCH_END_STATIC ON
      LINK_FLAGS ${client_static_link_flags})
  else()
    set_target_properties(${client_bin_name}
      PROPERTIES
      LINKER_LANGUAGE CXX
      LINK_SEARCH_END_STATIC ON)
  endif()
endfunction()
