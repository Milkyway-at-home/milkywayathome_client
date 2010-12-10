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

if(CMAKE_COMPILER_IS_GNUCC OR C_COMPILER_IS_CLANG)
  set(SSE_SSE2_FLAGS  "-msse2")
  set(SSE_SSE3_FLAGS  "-msse3 -msse2")
  set(SSE_SSSE3_FLAGS "-mssse3 -msse3 -msse2")
  set(SSE_SSE4_FLAGS  "-msse4 -mssse3 -msse3 -msse2")
  set(SSE_SSE41_FLAGS "-msse4.1 -msse4 -mssse3 -msse3 -msse2")
  set(SSE_SSE42_FLAGS "-msse4.2 -msse4.1 -msse4 -mssse3 -msse3 -msse2")
  set(SSE_SSE4A_FLAGS "-msse4a")

  set(SSE_LEVELS_X86_64 "sse3" "ssse3" "sse4" "sse41" "sse42")
  set(SSE_LEVELS_X86_32 "sse2")
  list(APPEND SSE_LEVELS_X86_32 ${SSE_LEVELS_X86_64})
elseif(MSVC)
  # Seems to not support anything beyond SSE2
  set(SSE_SSE2_FLAGS  "/arch:SSE2")
  set(SSE_LEVELS_X86_64 )
  set(SSE_LEVELS_X86_32 "sse2")
endif()

# x86_64 implies SSE and SSE2, so don't bother with special builds
if(SYSTEM_IS_X86_32)
  set(AVAILABLE_SSE_LEVELS ${SSE_LEVELS_X86_32} CACHE INTERNAL "Available SSE levels")
elseif(SYSTEM_IS_X86_64)
  set(AVAILABLE_SSE_LEVELS ${SSE_LEVELS_X86_64} CACHE INTERNAL "Available SSE levels")
endif()

function(add_sseX_library basename sse_level src headers)
  set(sse_lib_name ${basename}__${sse_level})
  add_library(${sse_lib_name} STATIC ${src} ${headers})

  string(TOUPPER ${sse_level} sse_level)
  set_target_properties(${sse_lib_name}
                          PROPERTIES
                            COMPILE_FLAGS "${SSE_${sse_level}_FLAGS}")
  set(SSE_LIB_NAME ${sse_lib_name} PARENT_SCOPE)
endfunction()

function(add_sseX_executable exe_basename sse_level src)
  set(sse_exe_name ${exe_basename}__${sse_level})
  add_executable(${sse_exe_name} ${src} ${headers})
  set_target_properties(${sse_exe_name}
                          PROPERTIES
                            COMPILE_FLAGS "${SSE_${sse_level}_FLAGS}")
  set(SSE_EXE_NAME ${sse_exe_name} PARENT_SCOPE)
endfunction()


