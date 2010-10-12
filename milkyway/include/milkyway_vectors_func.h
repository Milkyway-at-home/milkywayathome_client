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

#ifndef _MILKYWAY_VECTORS_FUNC_H_
#define _MILKYWAY_VECTORS_FUNC_H_

#include "real.h"
#include "milkyway_math_functions.h"


/* FIXME: Does not belong */
#define sqr(x) ((x) * (x))

__attribute__((const, always_inline))
inline mwvector mw_addv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    v.z = a.z + b.z;
    return v;
}

__attribute__((const, always_inline))
inline mwvector mw_subv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    v.z = a.z - b.z;

    return v;
}

__attribute__((const, always_inline))
inline mwvector mw_mulv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x * b.x;
    v.y = a.y * b.y;
    v.z = a.z * b.z;

    return v;
}

__attribute__((const, always_inline))
inline mwvector mw_divv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x / b.x;
    v.y = a.y / b.y;
    v.z = a.z / b.z;

    return v;
}

__attribute__((const, always_inline))
inline real mw_dotv(mwvector a, mwvector b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
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

    return mw_sqrt(sqr(a.x) + sqr(a.y) + sqr(a.z));
}

__attribute__((const, always_inline))
inline real mw_sqrv(mwvector a)
{
    return sqr(a.x) + sqr(a.y) + sqr(a.z);
}

__attribute__((const, always_inline))
inline real mw_absv(mwvector a)
{
    return mw_sqrt(mw_sqrv(a));
}

__attribute__((const, always_inline))
inline mwvector mw_mulvs(real s, mwvector a)
{
    mwvector v;
    v.x =  s * a.x;
    v.y =  s * a.y;
    v.z =  s * a.z;
    return v;
}

__attribute__((const, always_inline))
inline mwvector mw_mulmv(const mwmatrix m, mwvector a)
{
    mwvector tmp;
    tmp.x = mw_dotv(m[0], a);
    tmp.y = mw_dotv(m[1], a);
    tmp.z = mw_dotv(m[2], a);
    return tmp;
}

#define mw_incsubv(v1, v2) { (v1).x -= (v2).x; (v1).y -= (v2).y; (v1).z -= (v2).z; }
#define mw_incaddv(v1, v2) { (v1).x += (v2).x; (v1).y += (v2).y; (v1).z += (v2).z; }
#define mw_incdivv(v1, v2) { (v1).x /= (v2).x; (v1).y /= (v2).y; (v1).z /= (v2).z; }
#define mw_incmulv(v1, v2) { (v1).x *= (v2).x; (v1).y *= (v2).y; (v1).z *= (v2).z; }
#define mw_incnegv(v1) { (v1).x = -(v1).x; (v1).y = -(v1).y; (v1).z = -(v1).z; }
#define mw_zerov(v) { (v).x = 0.0; (v).y = 0.0; (v).z = 0.0; (v).w = 0.0; }

#define mw_incdivs(v, s) { (v).x /= (s); (v).y /= (s); (v).z /= (s); }
#define mw_incmulvs(v, s) { (v).x *= (s); (v).y *= (s); (v).z *= (s); }

#undef sqr


#endif /* _MILKYWAY_VECTORS_FUNC_H_ */

