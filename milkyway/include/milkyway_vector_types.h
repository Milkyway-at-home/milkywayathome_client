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

#include "real.h"

#ifndef _MSC_VER

typedef struct MW_ALIGN(4 * sizeof(real))
{
    real x, y, z, w;
} mwvector;

#else

typedef struct
{
    real x, y, z, w;
} mwvector;

#endif /* _MSC_VER */


#define L(v) ((v).x)
#define B(v) ((v).y)
#define R(v) ((v).z)

#define X(v) ((v).x)
#define Y(v) ((v).y)
#define Z(v) ((v).z)
#define W(v) ((v).w)


#define mw_vec(x, y, z) { (x), (y), (z), 0.0 }

#define SET_VECTOR(v, x, y, z) { X(v) = (x); Y(v) = (y); Z(v) = (z); }


#define NDIM 3

#define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
typedef mwvector mwmatrix[NDIM];
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

#define MWVECTOR_TYPE "Vector"
#define EMPTY_MWVECTOR { NAN, NAN, NAN, NAN }

#endif /* _MILKYWAY_VECTOR_TYPES_H_ */

