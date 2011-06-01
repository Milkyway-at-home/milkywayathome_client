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


macro(unset_cmake_default_dynamic)
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
endmacro()

if(NOT MSVC AND (CMAKE_BUILD_TYPE STREQUAL Release) AND NOT APPLE)
  set(strip_exe -s)
endif()

macro(maybe_static use_static)
  if(${use_static} MATCHES "ON")
    unset_cmake_default_dynamic()
    if(NOT MSVC)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread -static-libgcc -static-libstdc++")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -static-libgcc -static-libstdc++")
    endif()

  endif()
  if(NOT MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endmacro()


function(correct_static_link client_bin_name)
  if(UNIX AND NOT APPLE)
    set(client_static_link_flags "-static -static-libgcc -static-libstdc++")
  elseif(MINGW)
    set(client_static_link_flags "-static-libgcc -static-libstdc++")
  elseif(UNIX AND APPLE) # OS X
    #  No static
    set(client_static_link_flags "-static-libgcc -static-libstdc++")
  endif()


  set_target_properties(${client_bin_name}
                          PROPERTIES
                            LINKER_LANGUAGE CXX
                            LINK_FLAGS "${client_static_link_flags} ${strip_exe}"
                            LINK_SEARCH_END_STATIC ON)
endfunction()

function(milkyway_link client_bin_name use_boinc use_static link_libs)
  if(use_static)
    correct_static_link(${client_bin_name})
  else()
    if(use_boinc)
      if(NOT WIN32)
        list(APPEND link_libs "stdc++")
      endif()
    endif()

    set_target_properties(${client_bin_name}
                           PROPERTIES
                             LINKER_LANGUAGE CXX
                             LINK_FLAGS "${strip_exe}")
  endif()

  target_link_libraries(${client_bin_name} ${link_libs})
endfunction()

