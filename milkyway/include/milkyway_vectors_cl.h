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

#ifndef _MILKYWAY_VECTORS_CL_H_
#define _MILKYWAY_VECTORS_CL_H_

#include "milkyway_cl.h"
#include "real.h"
#include "milkyway_math_functions.h"

/* Wrap the basic operations on the kernel side */

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
inline real mw_absv(mwvector a)
{
    return mw_sqrt(mw_sqrv(a));
}

__attribute__((const, always_inline))
inline mwvector mw_mulvs(real s, mwvector a)
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
    mwvector tmp = mw_vec(mw_dotv(m[0], a),
                          mw_dotv(m[1], a),
                          mw_dotv(m[2], a));
    return tmp;
}


#define mw_incsubv(v1, v2) ((v1) -= (v2))
#define mw_incaddv(v1, v2) ((v1) += (v2))
#define mw_incdivv(v1, v2) ((v1) /= (v2))
#define mw_incmulv(v1, v2) ((v1) *= (v2))
#define mw_incnegv(v) ((v) = -(v))
#define mw_zerov(v) { (v).x = 0.0; (v).y = 0.0; (v).z = 0.0; (v).w = 0.0; }
#define mw_incdivs(v, s) { (v).x /= (s); (v).y /= (s); (v).z /= (s); (v).w /= (s) }
#define mw_incmulvs(v, s) { (v).x *= (s); (v).y *= (s); (v).z *= (s); (v).w *= (s) }

#endif /* _MILKYWAY_VECTORS_CL_H_ */

