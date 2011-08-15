# Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
# Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
# Magdon-Ismail and Rensselaer Polytechnic Institute.
#
# This file is part of Milkway@Home.
#
# Milkyway@Home is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Milkyway@Home is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
#

# CMAKE_SYSTEM_PROCESSOR is unfortunately still i386 on 64 bit OS X
if(    ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386"
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i486"
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i586"  # This is stupid
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i686"
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86"
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64"
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "amd64")
  set(SYSTEM_IS_X86 TRUE CACHE INTERNAL "Is x86")
else()
  set(SYSTEM_IS_X86 FALSE CACHE INTERNAL "Is x86")
endif()


if(   ${CMAKE_SYSTEM_PROCESSOR} MATCHES "powerpc"
    OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "ppc")
  set(SYSTEM_IS_PPC TRUE CACHE INTERNAL "Is PPC")
else()
  set(SYSTEM_IS_PPC FALSE CACHE INTERNAL "Is PPC")
endif()

# FIXME: This seems to not be set on windows
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(SYSTEM_IS_64 TRUE CACHE INTERNAL "Is 64 bit")
  set(SYSTEM_IS_32 FALSE CACHE INTERNAL "Is 32 bit")
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(SYSTEM_IS_64 FALSE CACHE INTERNAL "Is 64 bit")
  set(SYSTEM_IS_32 TRUE CACHE INTERNAL "Is 32 bit")
else()
  message(FATAL_ERROR "sizeof(void*) != 4, 8. What is this crazy system?")
endif()

if(SYSTEM_IS_X86 AND SYSTEM_IS_64)
  set(SYSTEM_IS_X86_64 TRUE CACHE INTERNAL "Is 64-bit x86")
endif()

if(SYSTEM_IS_X86 AND SYSTEM_IS_32)
  set(SYSTEM_IS_X86_32 TRUE CACHE INTERNAL "Is 32-bit x86")
endif()





