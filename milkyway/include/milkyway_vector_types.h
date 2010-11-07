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

#ifndef _MILKYWAY_VECTOR_TYPES_H_
#define _MILKYWAY_VECTOR_TYPES_H_

#include "milkyway_cl.h"
#include "real.h"

#if __OPENCL_VERSION__ >= __CL_VERSION_1_1__
  #define HAVE_VEC3 1
  #define USE_VEC3 1
#endif /* OpenCL 1.1 */

#if defined(__OPENCL_VERSION__)
 /* In the kernel */
  #if DOUBLEPREC
    typedef double2 real2;
    typedef double4 real4;
    #if HAVE_VEC3
    typedef double3 real3;
    #endif
  #else
    typedef float2 real2;
    typedef float4 real4;
  #if HAVE_VEC3
    typedef float3 real3;
  #endif
  #endif /* DOUBLEPREC */

  #if USE_VEC3
    typedef real3 mwvector;
  #else
    typedef real4 mwvector;
  #endif /* USE_VEC3 */
#else
  #if !defined(_MSC_VER) && !defined(__clang__)
/* When we align the structs we use the mwvector in, there is a clang
 * bug so for now we can't align it or use the ext_vector_type:
 * http://llvm.org/bugs/show_bug.cgi?id=8413 */
    #define MW_ALIGN(x) __attribute__((packed, aligned(x)))
  #else
    #define MW_ALIGN(x)
  #endif

typedef struct MW_ALIGN(4 * sizeof(real))
{
    real x, y, z, w;
} mwvector;


#endif /* __OPENCL_VERSION__ */


#define L(v) ((v).x)
#define B(v) ((v).y)
#define R(v) ((v).z)

#define X(v) ((v).x)
#define Y(v) ((v).y)
#define Z(v) ((v).z)
#define W(v) ((v).w)


#if USE_VEC3
  #define mw_vec(x, y, z) { (x), (y), (z) }
#else

/* Another glorious victory for standards */
#ifdef __NVIDIA_CL__
  #define mw_vec(x, y, z) ( (x), (y), (z), 0.0 )
#else
  #define mw_vec(x, y, z) { (x), (y), (z), 0.0 }
#endif /* __NVIDIA_CL__ */

#endif /* USE_VEC3 */


#define SET_VECTOR(v, x, y, z) { X(v) = (x); Y(v) = (y); Z(v) = (z); }



#define NDIM 3

#define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
typedef mwvector mwmatrix[NDIM];
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }


#endif /* _MILKYWAY_VECTOR_TYPES_H_ */

