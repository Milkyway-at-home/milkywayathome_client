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

include(CPUNameTest)

macro(set_os_specific_libs cl_required)

  get_info_from_processor_name()

  if(APPLE)
    if(${cl_required} MATCHES "ON")
      set(CMAKE_OSX_DEPLOYMENT_TARGET 10.6)
      set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.6.sdk")
    else()
      # Assume we only care about old OS X on PPC

      if(SYSTEM_IS_PPC)
        set(CMAKE_OSX_DEPLOYMENT_TARGET 10.3)
        set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.3.9.sdk")

        # You seem to have to specify these explicitly on really old OS X
        # CoreFoundation seems to take care of everything now
        find_library(CARBON_LIBRARY Carbon)
        find_library(STD_C_LIBRARY c)
        find_library(SYSTEM_STUBS SystemStubs)
        mark_as_advanced(CARBON_LIBRARY
                         SYSTEM_STUBS
                         STD_C_LIBRARY)
        list(APPEND OS_SPECIFIC_LIBS ${CARBON_LIBRARY} ${SYSTEM_STUBS} ${STD_C_LIBRARY})
      else()
        # Try to avoid the dyld: unknown required load command 0x80000022
        # runtime error on Leopard for binaries built on 10.6
        set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
        set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk")

        find_library(COREFOUNDATION_LIBRARY CoreFoundation)
        list(APPEND OS_SPECIFIC_LIBS ${COREFOUNDATION_LIBRARY})
      endif()
    endif()
endif()
  if(WIN32)
    set(OS_SPECIFIC_LIBS msvcrt)
  endif()
endmacro()

