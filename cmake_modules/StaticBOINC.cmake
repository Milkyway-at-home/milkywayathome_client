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

function(correct_static_link client_bin_name partially_dynamic)
  if(NOT APPLE)
    #CHECKME: What about Windows?
    set(client_static_link_flags "-static-libgcc -static-libstdc++")

    if(NOT partially_dynamic)
      # For dynamically linking to libOpenCL.so and dependencies
      set(client_static_link_flags "-static ${client_static_link_flags}")
    endif()

    if(UNIX AND NOT partially_dynamic)
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
    if(CMAKE_OSX_DEPLOYMENT_TARGET MATCHES "10.3")
      set(osx_link_flags "-bind_at_load")
    endif()

    set_target_properties(${client_bin_name}
                            PROPERTIES
                              LINKER_LANGUAGE CXX
                              LINK_FLAGS "${osx_link_flags}"
                              LINK_SEARCH_END_STATIC ON)
  endif()
endfunction()

function(milkyway_link client_bin_name use_boinc use_opencl use_static link_libs)
  if (UNIX AND NOT APPLE AND use_opencl)
    #Avoid weird linking conflict when trying to dynamically link
    #against opencl library, and avoid dynamically linking libc unless
    #we need to
    list(INSERT link_libs 0 c)
  endif()

  if(use_static)
    #Static linking tends to interfere with debugging with
    #valgrind/gdb in annoying ways

    # On Linux, you must dynamically link against libOpenCL.so, and
    #then pthreads, libc etc. to avoid conflicts. On Windows, you
    #statically link against libOpenCL.lib, which then links to a dll
    #at runtime.
    correct_static_link(${client_bin_name} ${use_opencl})
  else()
    if(use_boinc)
      list(APPEND link_libs "stdc++")
      if(UNIX)
        list(APPEND link_libs "pthread")
      endif()
    endif()
  endif()

  if(use_boinc)
    #Something dumb is happening with linking order and I'm sick of
    #fighting with it right now (Sep. 2010)
    list(APPEND link_libs ${BOINC_LIBRARIES})
    list(INSERT link_libs 0 ${BOINC_LIBRARIES})
  endif()

  #Include libm first before anything else. I can't seem to make cmake
  #always give a static executable, and only try to use static
  #libraries. When statically linking, by explicitly finding libm.a,
  #and linking to that first, it avoids trying to dynamically link
  #libm

  target_link_libraries(${client_bin_name} ${link_libs})
endfunction()

macro(maybe_static)
  if(STATIC_EVERYTHING)
    if(NOT APPLE)
      set(LIBM_USE_STATIC 1)
    endif()

    set(POPT_USE_STATIC 1)
    set(JSON_C_USE_STATIC 1)
    set(LIBINTL_USE_STATIC 1)
    set(LIBICONV_USE_STATIC 1)
    set(BOINC_USE_STATIC 1)
  endif()
endmacro()


