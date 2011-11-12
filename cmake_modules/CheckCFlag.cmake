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

include(CheckCCompilerFlag)

function(add_flag_if_supported flagname)
  check_c_compiler_flag("${flagname}" HAVE_FLAG_${flagname})

  if(${HAVE_FLAG_${flagname}})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flagname}" PARENT_SCOPE)
  endif()
endfunction()

function(append_supported_flags flags)
  foreach(flag ${flags})
    add_flag_if_supported(${flag})
  endforeach()
endfunction()


check_c_compiler_flag("/fp:fast" HAVE_FP_FAST_MATH)
if(HAVE_FP_FAST_MATH)
  set(FAST_MATH_FLAGS "/fp:fast")
endif()

check_c_compiler_flag("-ffast-math" HAVE_F_FAST_MATH)
if(HAVE_F_FAST_MATH)
  set(FAST_MATH_FLAGS "-ffast-math")
endif()


#check_c_compiler_flag("-march=em64t" HAVE_MARCH_EM64T)
#check_c_compiler_flag("-march=anyx86" HAVE_MARCH_ANYX86)
#check_c_compiler_flag("-mtune=generic" HAVE_MTUNE_GENERIC)

check_c_compiler_flag("-static-libstdc++" HAVE_FLAG_STATIC_LIBSTDCPP)
check_c_compiler_flag("-static-libgcc" HAVE_FLAG_STATIC_LIBGCC)
check_c_compiler_flag("-pthread" HAVE_FLAG_PTHREAD)

