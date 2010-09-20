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

macro(unknown_system)
  message(WARNING "Unknown system: ${CMAKE_SYSTEM}")
  set(boinc_sys_name "unknown-system")
endmacro()

function(get_boinc_bin_name basename version plan)
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")
    set(SYSTEM_IS_X86 TRUE)
  else()
    set(SYSTEM_IS_X86 FALSE)
  endif()

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ppc")
    set(SYSTEM_IS_PPC TRUE)
  else()
    set(SYSTEM_IS_PPC FALSE)
  endif()

  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(SYSTEM_IS_64 TRUE)
  elseif(PTR_SIZE EQUAL 4)
    set(SYSTEM_IS_64 FALSE)
  else()
    message(FATAL_ERROR "sizeof(void*) != 4, 8. What is this crazy system?")
  endif()

  if(WIN32)
    if(NOT SYSTEM_IS_X86)
      message("Weird Windows that only sort of once existed?")
    endif()

    if(SYSTEM_IS_64)
      set(boinc_sys_name "windows_intelx86")
    else()
      set(boinc_sys_name "windows_x86_64")
    endif()
  elseif(UNIX)
    if(APPLE)
      if(SYSTEM_IS_64 AND SYSTEM_IS_X86)
        set(boinc_sys_name "x86_64-apple-darwin")
      elseif(NOT SYSTEM_IS_64 AND SYSTEM_IS_x86)
        set(boinc_sys_name "i686-apple-darwin")
      elseif(SYSTEM_IS_PPC)
        set(boinc_sys_name "powerpc-apple-darwin")
      else()
        unknown_system()
      endif()
    endif()

    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      if(SYSTEM_IS_X86 AND SYSTEM_IS_64)
        set(boinc_sys_name "x86_64-pc-linux-gnu")
      elseif(SYSTEM_IS_X86 AND NOT SYSTEM_IS_64)
        set(boinc_sys_name "i686-pc-linux-gnu")
      elseif(SYSTEM_IS_PPC AND SYSTEM_IS_64)
        set(boinc_sys_name "ppc64-linux-gnu")
      elseif(SYSTEM_IS_PPC AND NOT SYSTEM_IS_64)
        set(boinc_sys_name "powerpc-linux-gnu")
      else()
        unknown_system()
      endif()
    endif()

    if(${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
      if(SYSTEM_IS_X86 AND SYSTEM_IS_64)
        set(boinc_sys_name "x86_64-pc-freebsd")
      elseif(SYSTEM_IS_X86 AND NOT SYSTEM_IS_64)
        set(boinc_sys_name "i686-pc-freebsd")
      else()
        unknown_system()
      endif()
    endif()
  else()
    unknown_system()
  endif()

  if(plan)
    if(plan STREQUAL "sse2" AND SYSTEM_IS_64 AND SYSTEM_IS_X86)
      # x86_64 implies sse2
      set(boinc_suffix "")
    else()
      set(boinc_suffix "__${plan}")
    endif()
  else()
    # No __ if no plan set
    set(boinc_suffix "")
  endif()

  set(BOINC_BIN_NAME "${basename}_${version}_${boinc_sys_name}${boinc_suffix}" PARENT_SCOPE)
endfunction()


