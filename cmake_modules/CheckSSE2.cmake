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
include(CheckIncludeFiles)
include(CheckCCompilerFlag)

macro(str_append xs x)
  set(${xs} "${${xs}} ${x}")
endmacro()


if(SYSTEM_IS_X86)
  if(MSVC OR "${CMAKE_C_COMPILER_ID}" MATCHES "SunPro")
    set(NEED_SSE_DEFINES TRUE)
  endif()

  if("${CMAKE_C_COMPILER_ID}" MATCHES "SunPro")
    set(HAVE_FLAG_XARCH_SSE TRUE)
    set(HAVE_FLAG_XARCH_SSE2 TRUE)
    set(HAVE_FLAG_XARCH_SSE3 TRUE)
    set(HAVE_FLAG_XARCH_SSE4 TRUE)
    set(HAVE_FLAG_XARCH_SSE4_1 TRUE)
  endif()

  if(NOT MSVC)
    check_c_compiler_flag("-mfpmath=sse" HAVE_FLAG_M_FPMATH_SSE)
    check_c_compiler_flag("-msse" HAVE_FLAG_M_SSE)
    check_c_compiler_flag("-msse2" HAVE_FLAG_M_SSE2)
    check_c_compiler_flag("-msse3" HAVE_FLAG_M_SSE3)
    check_c_compiler_flag("-msse4" HAVE_FLAG_M_SSE4)
    check_c_compiler_flag("-msse4.1" HAVE_FLAG_M_SSE41)
    check_c_compiler_flag("-mavx" HAVE_FLAG_M_AVX)


    # These all fail for some reason
    #check_c_compiler_flag("-xarch=sse" HAVE_FLAG_XARCH_SSE)
    #check_c_compiler_flag("-xarch=sse2" HAVE_FLAG_XARCH_SSE2)
    #check_c_compiler_flag("-xarch=sse3" HAVE_FLAG_XARCH_SSE3)
    #check_c_compiler_flag("-xarch=sse4" HAVE_FLAG_XARCH_SSE4)
    #check_c_compiler_flag("-xarch=sse4_1" HAVE_FLAG_XARCH_SSE4_1)


    if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
      # Check succeeds when it shouldn't for some reason.
      set(HAVE_FLAG_M_SSE4_1 FALSE)
    else()
      check_c_compiler_flag("-msse4_1" HAVE_FLAG_M_SSE4_1) # pathcc uses sse4_1 names
    endif()

    set(BASE_SSE_FLAGS )
    if(HAVE_FLAG_M_FPMATH_SSE)
      str_append(BASE_SSE_FLAGS "-mfpmath=sse")
    endif()
    if(HAVE_FLAG_M_SSE)
      str_append(BASE_SSE_FLAGS "-msse")
    endif()
    if(HAVE_FLAG_XARCH_SSE)
      str_append(BASE_SSE_FLAGS "-xarch=sse")
    endif()

    set(SSE2_FLAGS ${BASE_SSE_FLAGS})
    if(HAVE_FLAG_M_SSE2)
      str_append(SSE2_FLAGS "-msse2")
    endif()
    if(HAVE_FLAG_XARCH_SSE2)
      str_append(SSE2_FLAGS "-xarch=sse2")
    endif()

    set(SSE3_FLAGS ${SSE2_FLAGS})
    if(HAVE_FLAG_M_SSE3)
      str_append(SSE3_FLAGS "-msse3")
    endif()
    if(HAVE_FLAG_XARCH_SSE3)
      str_append(SSE3_FLAGS "-xarch=sse3")
    endif()

    set(SSE4_FLAGS ${SSE3_FLAGS})
    if(HAVE_FLAG_M_SSE4)
      str_append(SSE4_FLAGS "-msse4")
    endif()
    if(HAVE_FLAG_XARCH_SSE4)
      str_append(SSE4_FLAGS "-xarch=sse4")
    endif()

    set(SSE41_FLAGS ${SSE4_FLAGS})
    if(HAVE_FLAG_M_SSE41)
      str_append(SSE41_FLAGS "-msse4.1")
    endif()
    if(HAVE_FLAG_M_SSE4_1)
      str_append(SSE41_FLAGS "-msse4_1")
    endif()
    if(HAVE_FLAG_XARCH_SSE4_1)
      str_append(SSE41_FLAGS "-xarch=sse4_1")
    endif()


    set(AVX_FLAGS ${SSE41_FLAGS})
    if(HAVE_FLAG_M_AVX)
      str_append(AVX_FLAGS "-mavx")
    endif()
    if(HAVE_FLAG_XARCH_AVX)
      str_append(AVX_FLAGS "-xarch=avx")
    endif()


    check_c_compiler_flag("-mfpmath=387" HAVE_FLAG_M_FPMATH_387)
    check_c_compiler_flag("-mno-sse" HAVE_FLAG_M_NO_SSE)
    check_c_compiler_flag("-mno-sse2" HAVE_FLAG_M_NO_SSE2)
    check_c_compiler_flag("-mno-sse3" HAVE_FLAG_M_NO_SSE3)
    check_c_compiler_flag("-mno-sse4" HAVE_FLAG_M_NO_SSE4)
    check_c_compiler_flag("-mno-sse4.1" HAVE_FLAG_M_NO_SSE41)
    check_c_compiler_flag("-mno-sse4_1" HAVE_FLAG_M_NO_SSE4_1)
    check_c_compiler_flag("-mno-avx" HAVE_FLAG_M_NO_AVX)

    set(DISABLE_SSE_FLAGS )
    if(HAVE_FLAG_M_FPMATH_387)
      str_append(DISABLE_SSE_FLAGS "-mfpmath=387")
    endif()
    if(HAVE_FLAG_M_NO_SSE)
      str_append(DISABLE_SSE_FLAGS "-mno-sse")
    endif()

    set(DISABLE_SSE2_FLAGS ${DISABLE_SSE_FLAGS})
    if(HAVE_FLAG_M_NO_SSE2)
      str_append(DISABLE_SSE2_FLAGS "-mno-sse2")
    endif()

    set(DISABLE_SSE3_FLAGS "")
    if(HAVE_FLAG_M_NO_SSE3)
      str_append(DISABLE_SSE3_FLAGS "-mno-sse3")
    endif()

    set(DISABLE_SSE4_FLAGS "")
    if(HAVE_FLAG_M_NO_SSE4)
      str_append(DISABLE_SSE4_FLAGS "-mno-sse4")
    endif()

    set(DISABLE_SSE41_FLAGS "")
    if(HAVE_FLAG_M_NO_SSE41)
      str_append(DISABLE_SSE41_FLAGS "-mno-sse4.1")
    endif()
    if(HAVE_FLAG_M_NO_SSE4_1)
      str_append(DISABLE_SSE41_FLAGS "-mno-sse4_1")
    endif()

    set(DISABLE_AVX_FLAGS "")
    if(HAVE_FLAG_M_NO_AVX)
      str_append(DISABLE_AVX_FLAGS "-mno-avx")
    endif()
  else()
    #set(BASE_SSE_FLAGS "/arch:SSE")
    set(SSE2_FLAGS "/arch:SSE2")
    set(DISABLE_SSE2_FLAGS "")
    set(DISABLE_SSE3_FLAGS "")
    set(DISABLE_SSE41_FLAGS "")
    set(DISABLE_AVX_FLAGS "")

    # MSVC doesn't generate SSE3 itself, and doesn't define this
    set(SSE3_FLAGS "${SSE2_FLAGS}")
    set(SSE41_FLAGS "${SSE3_FLAGS}")
    set(AVX_FLAGS "/arch:AVX")
  endif()

  if(NEED_SSE_DEFINES)
    str_append(BASE_SSE_FLAGS "-D__SSE__=1")
    str_append(SSE2_FLAGS "-D__SSE2__=1")
    str_append(SSE3_FLAGS "-D__SSE3__=1")
    str_append(SSE4_FLAGS "-D__SSE4__=1")
    str_append(SSE41_FLAGS "-D__SSE4_1__=1")
    str_append(AVX_FLAGS "-D__AVX__=1")

    # want the defines but not /arch:SSE2
    str_append(AVX_FLAGS "-D__SSE4_1__=1")
    str_append(AVX_FLAGS "-D__SSE3__=1")
    str_append(AVX_FLAGS "-D__SSE2__=1")
  endif()
endif()


# On OS X 10.6 the macports GCC can support AVX, but you must use the
# system assembler which doesn't and fails
set(_CMAKE_C_FLAGS "${CMAKE_C_FLAGS}") # Other methods of passing cflags to try_compile don't seem to work
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${AVX_FLAGS}")
try_compile(AVX_CHECK ${CMAKE_BINARY_DIR} ${MILKYWAYATHOME_CLIENT_CMAKE_MODULES}/test_avx.c)
set(CMAKE_C_FLAGS ${_CMAKE_C_FLAGS})
if(AVX_CHECK)
  message(STATUS "AVX compiler flags - '${AVX_FLAGS}'")
  set(HAVE_AVX TRUE CACHE INTERNAL "Compiler has AVX support")
endif()
mark_as_advanced(HAVE_AVX)


set(CMAKE_REQUIRED_FLAGS "${SSE41_FLAGS}")
check_include_files(smmintrin.h HAVE_SSE41 CACHE INTERNAL "Compiler has SSE4.1 headers")
mark_as_advanced(HAVE_SSE41)

set(CMAKE_REQUIRED_FLAGS "${SSE3_FLAGS}")
check_include_files(pmmintrin.h HAVE_SSE3 CACHE INTERNAL "Compiler has SSE3 headers")
mark_as_advanced(HAVE_SSE3)

set(CMAKE_REQUIRED_FLAGS "${SSE2_FLAGS}")
check_include_files(emmintrin.h HAVE_SSE2 CACHE INTERNAL "Compiler has SSE2 headers")
mark_as_advanced(HAVE_SSE2)

set(CMAKE_REQUIRED_FLAGS "")

if(APPLE AND SYSTEM_IS_X86)
  set(ALWAYS_HAVE_SSE2 TRUE CACHE INTERNAL "System always has SSE2")
elseif(SYSTEM_IS_X86_64)
  set(ALWAYS_HAVE_SSE2 TRUE CACHE INTERNAL "System always has SSE2")
else()
  set(ALWAYS_HAVE_SSE2 FALSE CACHE INTERNAL "System always has SSE2")
endif()

if(APPLE AND SYSTEM_IS_X86)
  set(ALWAYS_HAVE_SSE3 TRUE CACHE INTERNAL "System always has SSE3")
else()
  set(ALWAYS_HAVE_SSE3 FALSE CACHE INTERNAL "System always has SSE3")
endif()

function(disable_sse41 target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()
  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${DISABLE_SSE41_FLAGS}")
endfunction()

function(disable_sse3 target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()
  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${DISABLE_SSE3_FLAGS}")
endfunction()

function(disable_sse2 target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()
  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${DISABLE_SSE2_FLAGS}")
endfunction()

function(enable_sse41 target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()

  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${SSE41_FLAGS}")
  get_target_property(new_comp_flags ${target} COMPILE_FLAGS)
endfunction()

function(enable_sse3 target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()

  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${SSE3_FLAGS}")
  get_target_property(new_comp_flags ${target} COMPILE_FLAGS)
endfunction()

function(enable_sse2 target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()

  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${SSE2_FLAGS}")
endfunction()

function(enable_avx target)
  get_target_property(comp_flags ${target} COMPILE_FLAGS)
  if(comp_flags STREQUAL "comp_flags-NOTFOUND")
    set(comp_flags "")
  endif()

  set_target_properties(${target}
                          PROPERTIES
                            COMPILE_FLAGS "${comp_flags} ${AVX_FLAGS}")
endfunction()


function(maybe_disable_ssen)
  if(SYSTEM_IS_X86)
    foreach(i ${ARGN})
      if(NOT ALWAYS_HAVE_SSE2)
        disable_sse2(${i})
      endif()
      if(NOT ALWAYS_HAVE_SSE3)
        disable_sse3(${i})
      endif()
    endforeach()
  endif()
endfunction()

