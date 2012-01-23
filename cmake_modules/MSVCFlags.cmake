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

include(CPUNameTest)

macro(set_msvc_extra_link_flags)
  if(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "/VERBOSE:LIB ${CMAKE_EXE_LINKER_FLAGS}")

    # Minimum is 5.0 on x86 and (5.02 on x64, same as default)
    if(SYSTEM_IS_32)
      set(CMAKE_EXE_LINKER_FLAGS "/SUBSYSTEM:CONSOLE,5.0 ${CMAKE_EXE_LINKER_FLAGS}")
    endif()

    if(MSVC_USE_STATIC_CRT)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:msvcrt.lib")
    endif()
    # Do not use /GL (whole program optimization)
    set(_extra_msvc_opt_flags "/fp:precise /fp:except- /Ot /Oi /Ob2 ${CMAKE_C_FLAGS}")
    set(CMAKE_C_FLAGS "${_extra_msvc_opt_flags} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${_extra_msvc_opt_flags} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_C_FLAGS_RELEASE "${_extra_msvc_opt_flags} ${CMAKE_C_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELEASE "${_extra_msvc_opt_flags} ${CMAKE_CXX_FLAGS_RELEASE}")
  endif()
endmacro()


macro(set_msvc_mt)
  foreach(flag_var
      CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO

      CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
      CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
    if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    endif(${flag_var} MATCHES "/MD")
  endforeach(flag_var)
endmacro()


macro(set_msvc_find_paths)
  set(MSVC_FIND_ROOT "${CMAKE_BINARY_DIR}/../MSVCDeps")

  if(MSVC_USE_STATIC_CRT)
    set(MT_OR_MD "mt")
  else()
    set(MT_OR_MD "md")
  endif()

  set(MSVC_ARCH "Win32")
  set(MSVC_BUILD "Release")

  list(APPEND CMAKE_LIBRARY_PATH ${CMAKE_BINARY_DIR}/lib/${MSVC_BUILD})
  list(APPEND CMAKE_INCLUDE_PATH ${CMAKE_BINARY_DIR}/include)

  list(APPEND CMAKE_LIBRARY_PATH "${MSVC_FIND_ROOT}/lib/${MSVC_BUILD}/${MSVC_ARCH}/${MT_OR_MD}")
  list(APPEND CMAKE_INCLUDE_PATH "${MSVC_FIND_ROOT}/include")
endmacro()


