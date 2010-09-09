/*
Copyright 2010 Matthew Arsenault

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

#ifndef _MILKYWAY_VECTORS_CLANG_H_
#define _MILKYWAY_VECTORS_CLANG_H_

#include "real.h"
#include "milkyway_math_functions.h"

/* Clang's vector stuff is nicer and what's in OpenCL */

typedef real mwvector __attribute__((ext_vector_type(4)));

#define L(v) ((v).x)
#define B(v) ((v).y)
#define R(v) ((v).z)

#define X(v) ((v).x)
#define Y(v) ((v).y)
#define Z(v) ((v).z)


#define mw_vec(x, y, z) { (x), (y), (z), 0.0 }

#define NDIM 3

typedef mwvector mwmatrix[NDIM];

#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }


__attribute__((const, always_inline))
inline mwvector mw_addv(mwvector a, mwvector b)
{
    return a + b;
}

__attribute__((const, always_inline))
inline mwvector mw_subv(mwvector a, mwvector b)
{
    return a - b;
}

__attribute__((const, always_inline))
inline mwvector mw_mulv(mwvector a, mwvector b)
{
    return a * b;
}

__attribute__((const, always_inline))
inline mwvector mw_divv(mwvector a, mwvector b)
{
    return a / b;
}

__attribute__((const, always_inline))
inline real mw_dotv(mwvector a, mwvector b)
{
    mwvector tmp = a * b;
    return tmp.x + tmp.y + tmp.z + tmp.w;
}

__attribute__((const, always_inline))
inline mwvector mw_crossv(mwvector a, mwvector b)
{
    mwvector tmp = mw_vec(b.z * a.y - b.y * a.z,
                          b.x * a.z - b.z * a.x,
                          b.y * a.x - b.x * a.y);
    return tmp;
}

__attribute__((const, always_inline))
inline real mw_length(mwvector a)
{
    mwvector tmp = a * a;
    return mw_sqrt(tmp.x + tmp.y + tmp.z + tmp.w);
}

__attribute__((const, always_inline))
inline real mw_sqrv(mwvector a)
{
    return mw_dotv(a, a);
}

__attribute__((const, always_inline))
inline mwvector mw_mulvs(mwvector a, real s)
{
    mwvector tmp =
        {
            s * a.x,
            s * a.y,
            s * a.z,
            s * a.w
        };
    return tmp;
}

__attribute__((pure, const, always_inline))
inline mwvector mw_mulmv(const mwmatrix m, mwvector a)
{
    mwvector tmp = mw_vec(mw_dotv(m[0], a)
                          mw_dotv(m[1], a),
                          mw_dotv(m[2], a));
    return tmp;
}


#endif /* _MILKYWAY_VECTORS_CLANG_H_ */

