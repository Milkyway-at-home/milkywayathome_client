/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
  #error "Only milkyway_math.h can be included directly."
#endif

#ifndef _MILKYWAY_VECTORS_H_
#define _MILKYWAY_VECTORS_H_


#include "real.h"

#if NBODY_OPENCL || SEPARATION_OPENCL || defined(__OPENCL_VERSION__)
  #include "milkyway_vectors_cl.h"
#else
  #include "milkyway_vectors_cpu.h"
#endif

#include "milkyway_vectors_ops.h"


#if 0
  #if defined(__clang__) || defined(__OPENCL_VERSION__)
    #include "milkyway_vectors_clang.h"
  #else /* GCC, MSVC */
    #include "milkyway_vectors_gcc.h"
  #endif /* defined(__clang__) || defined(__OPENCL_VERSION__) */
#endif /* 0 */

#endif /* _MILKYWAY_VECTORS_H_ */

