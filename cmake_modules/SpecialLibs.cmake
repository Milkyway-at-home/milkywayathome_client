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

if(APPLE)
  if(IS_DIRECTORY "/Developer/SDKs/MacOSX10.3.9.sdk")
    set(HAVE_10_3_SDK TRUE)
  endif()
  # if(IS_DIRECTORY "/Developer/SDKs/MacOSX10.4u.sdk")
  #   set(HAVE_10_4_SDK TRUE)
  # endif()
  # if(IS_DIRECTORY "/Developer/SDKs/MacOSX10.5.sdk")
  #   set(HAVE_10_5_SDK TRUE)
  # endif()
  # if(IS_DIRECTORY "/Developer/SDKs/MacOSX10.6.sdk")
  #   set(HAVE_10_6_SDK TRUE)
  # endif()
  # if(IS_DIRECTORY "/Developer/SDKs/MacOSX10.7.sdk")
  #   set(HAVE_10_7_SDK TRUE)
  # endif()

  if(MILKYWAY_OPENCL)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.6 CACHE STRING "" FORCE)
#    set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.6.sdk" CACHE PATH "" FORCE)

    find_library(COREFOUNDATION_LIBRARY CoreFoundation)
    list(APPEND OS_SPECIFIC_LIBS ${COREFOUNDATION_LIBRARY})
  else()
    # Assume we only care about old OS X on PPC
    if(SYSTEM_IS_PPC)
      if (NOT HAVE_10_3_SDK)
        message(FATAL "OS X 10.3 SDK required for PPC build")
      endif()

      set(CMAKE_OSX_DEPLOYMENT_TARGET 10.3 CACHE STRING "" FORCE)
      set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.3.9.sdk" CACHE PATH "" FORCE)

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

      if(SYSTEM_IS_32)
        set(CMAKE_OSX_DEPLOYMENT_TARGET 10.4 CACHE STRING "" FORCE)
      else()
        set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5 CACHE STRING "" FORCE)
      endif()

      #set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.4u.sdk" CACHE PATH "" FORCE)

      # if(CMAKE_OSX_ARCHITECTURES MATCHES "i386")
      #   if(NOT HAVE_10_4_SDK)
      #     message(FATAL "OS X 10.4 SDK required for x86_32 build")
      #   endif()
      #   # When building for 32 bit, seem to need the 10.4 SDK
      #   set(CMAKE_OSX_DEPLOYMENT_TARGET 10.4 CACHE STRING "" FORCE)
      #   #set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.4u.sdk" CACHE PATH "" FORCE)
      # else()
      #   # Try for the lowest version SDK we can
      #   if(HAVE_10_5_SDK)
      #     set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5 CACHE STRING "" FORCE)
      #     set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk" CACHE PATH "" FORCE)
      #   elseif(HAVE_10_6_SDK)
      #     set(CMAKE_OSX_DEPLOYMENT_TARGET 10.6 CACHE STRING "" FORCE)
      #    #set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.6.sdk" CACHE PATH "" FORCE)
      #   elseif(HAVE_10_7_SDK)
      #     set(CMAKE_OSX_DEPLOYMENT_TARGET 10.7 CACHE STRING "" FORCE)
      #    #set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.7.sdk" CACHE PATH "" FORCE)
      #   else()
      #     message(FATAL "No OS X SDKs found")
      #   endif()
      # endif()
      find_library(COREFOUNDATION_LIBRARY CoreFoundation)
      list(APPEND OS_SPECIFIC_LIBS ${COREFOUNDATION_LIBRARY})
    endif()
  endif()

  if(NOT IS_DIRECTORY ${CMAKE_OSX_SYSROOT})
    message(FATAL_ERROR "Correct OS X SDK installation version missing: ${CMAKE_OSX_SYSROOT}")
  endif()
endif()

#if(WIN32)
#  if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
#    set(OS_SPECIFIC_LIBS msvcrtd)
#  else()
#    set(OS_SPECIFIC_LIBS msvcrt)
#  endif()
#endif()



