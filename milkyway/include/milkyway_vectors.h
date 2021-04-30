/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#if !defined(_MILKYWAY_MATH_H_INSIDE_) && !defined(MILKYWAY_MATH_COMPILATION)
  #error "Only milkyway_math.h can be included directly."
#endif

#ifndef _MILKYWAY_VECTORS_H_
#define _MILKYWAY_VECTORS_H_

typedef struct MW_ALIGN_TYPE_V(4 * sizeof(real))
{
    real x, y, z, w;
} mwvector;

#define L(v) ((v).x)
#define B(v) ((v).y)
#define R(v) ((v).z)

#define X(v) ((v).x)
#define Y(v) ((v).y)
#define Z(v) ((v).z)
#define W(v) ((v).w)

/* Cylinderical Macros */
#define CR(v) ((v).x)
#define CT(v) ((v).y)
#define CZ(v) ((v).z)


#define mw_vec(x, y, z) { (x), (y), (z), 0.0 }

#define SET_VECTOR(v, x, y, z) { X(v) = (x); Y(v) = (y); Z(v) = (z); }


#define NDIM 3

#define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
typedef mwvector mwmatrix[NDIM];
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

#define MWVECTOR_TYPE "Vector"


CONST_F
static inline mwvector mw_addv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    v.z = a.z + b.z;
    return v;
}

CONST_F
static inline mwvector mw_subv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    v.z = a.z - b.z;

    return v;
}

CONST_F
static inline mwvector mw_mulv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x * b.x;
    v.y = a.y * b.y;
    v.z = a.z * b.z;

    return v;
}

CONST_F
static inline mwvector mw_divv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = a.x / b.x;
    v.y = a.y / b.y;
    v.z = a.z / b.z;

    return v;
}

CONST_F
static inline real mw_dotv(mwvector a, mwvector b)
{
    return mw_mad(a.z, b.z, mw_mad(a.y, b.y, a.x * b.x));
}

CONST_F
static inline mwvector mw_crossv(mwvector a, mwvector b)
{
    mwvector tmp = mw_vec(b.z * a.y - b.y * a.z,
                          b.x * a.z - b.z * a.x,
                          b.y * a.x - b.x * a.y);
    return tmp;
}

CONST_F
static inline real mw_length(mwvector a)
{
    return mw_sqrt(sqr(a.x) + sqr(a.y) + sqr(a.z));
}

CONST_F
static inline real mw_sqrv(mwvector a)
{
    return mw_mad(a.z, a.z, mw_mad(a.y, a.y, a.x * a.x));
}

CONST_F
static inline real mw_absv(mwvector a)
{
    return mw_sqrt(mw_sqrv(a));
}

CONST_F
static inline mwvector mw_mulvs(mwvector a, real s)
{
    mwvector v;
    v.x =  s * a.x;
    v.y =  s * a.y;
    v.z =  s * a.z;
    return v;
}

CONST_F
static inline mwvector mw_divvs(mwvector a, real s)
{
    mwvector v;
    v.x = a.x / s;
    v.y = a.y / s;
    v.z = a.z / s;
    return v;
}

CONST_F
static inline mwvector mw_negv(mwvector a)
{
    mwvector v;
    v.x = -a.x;
    v.y = -a.y;
    v.z = -a.z;
    return v;
}

CONST_F
static inline mwvector mw_mulmv(const mwmatrix m, mwvector a)
{
    mwvector tmp;
    tmp.x = mw_dotv(m[0], a);
    tmp.y = mw_dotv(m[1], a);
    tmp.z = mw_dotv(m[2], a);
    return tmp;
}

CONST_F
static inline real mw_distv(mwvector u, mwvector v)
{
    return mw_sqrt(sqr(u.x - v.x) + sqr(u.y - v.y) + sqr(u.z - v.z));
}

/* Angle between two vectors, in the range [0,pi] */
CONST_F
static inline real mw_vecangle(mwvector a, mwvector b)
{
    real anorm, bnorm, dot;

    anorm = mw_length(a);
    bnorm = mw_length(b);
    dot = mw_dotv(a, b);

    return mw_acos(dot / (anorm * bnorm));
}


#define mw_incsubv(v1, v2) { (v1).x -= (v2).x; (v1).y -= (v2).y; (v1).z -= (v2).z; }
#define mw_incaddv(v1, v2) { (v1).x += (v2).x; (v1).y += (v2).y; (v1).z += (v2).z; }
#define mw_incdivv(v1, v2) { (v1).x /= (v2).x; (v1).y /= (v2).y; (v1).z /= (v2).z; }
#define mw_incmulv(v1, v2) { (v1).x *= (v2).x; (v1).y *= (v2).y; (v1).z *= (v2).z; }
#define mw_incnegv(v1) { (v1).x = -(v1).x; (v1).y = -(v1).y; (v1).z = -(v1).z; }
#define mw_zerov(v) { (v).x = 0.0; (v).y = 0.0; (v).z = 0.0; }

/* v1 += s * v2 */
#define mw_incaddv_s(v1, v2, s) { (v1).x += (s) * (v2).x; (v1).y += (s) * (v2).y; (v1).z += (s) * (v2).z; }

#if !USE_MAD

/* v1 -= s * v2 */
#define mw_incsubv_s(v1, v2, s) { (v1).x -= (s) * (v2).x; (v1).y -= (s) * (v2).y; (v1).z -= (s) * (v2).z; }

#else

#define mw_incsubv_s(v1, v2, s)                 \
    {                                           \
        (v1).x = mw_mad(-(s), (v2).x, (v1).x);  \
        (v1).y = mw_mad(-(s), (v2).y, (v1).y);  \
        (v1).z = mw_mad(-(s), (v2).z, (v1).z);  \
    }

#endif /* USE_MAD */

#define mw_incdivs(v, s) { (v).x /= (s); (v).y /= (s); (v).z /= (s); }
#define mw_incmulvs(v, s) { (v).x *= (s); (v).y *= (s); (v).z *= (s); }

#define mw_normalize(v) { real len = mw_length(v); (v).x /= len; (v).y /= len; (v).z /= len; }

/* Outer product */

static inline void mw_outv(mwmatrix p, mwvector v, mwvector u)
{
    X(p[0]) = X(v) * X(u);
    Y(p[0]) = X(v) * Y(u);
    Z(p[0]) = X(v) * Z(u);

    X(p[1]) = Y(v) * X(u);
    Y(p[1]) = Y(v) * Y(u);
    Z(p[1]) = Y(v) * Z(u);

    X(p[2]) = Z(v) * X(u);
    Y(p[2]) = Z(v) * Y(u);
    Z(p[2]) = Z(v) * Z(u);
}

/* Outer product of a vector with itself */
static inline void mw_outsqrv(mwmatrix p, mwvector v)
{
    X(p[0]) = X(v) * X(v);
    Y(p[0]) = X(v) * Y(v);
    Z(p[0]) = X(v) * Z(v);

    X(p[1]) = Y(p[0]);
    Y(p[1]) = Y(v) * Y(v);
    Z(p[1]) = Y(v) * Z(v);

    X(p[2]) = Z(p[0]);
    Y(p[2]) = Z(p[1]);
    Z(p[2]) = Z(v) * Z(v);
}

static inline void mw_set_diagonal_matrix(mwmatrix p, real s)
{
    X(p[0]) = s;
    Y(p[0]) = 0.0;
    Z(p[0]) = 0.0;

    X(p[1]) = 0.0;
    Y(p[1]) = s;
    Z(p[1]) = 0.0;

    X(p[2]) = 0.0;
    Y(p[2]) = 0.0;
    Z(p[2]) = s;
}

static inline void mw_set_matrix_identity(mwmatrix p)
{
    X(p[0]) = 1.0;
    Y(p[0]) = 0.0;
    Z(p[0]) = 0.0;

    X(p[1]) = 0.0;
    Y(p[1]) = 1.0;
    Z(p[1]) = 0.0;

    X(p[2]) = 0.0;
    Y(p[2]) = 0.0;
    Z(p[2]) = 1.0;
}

static inline void mw_set_matrix_zero(mwmatrix p)
{
    X(p[0]) = 0.0;
    Y(p[0]) = 0.0;
    Z(p[0]) = 0.0;

    X(p[1]) = 0.0;
    Y(p[1]) = 0.0;
    Z(p[1]) = 0.0;

    X(p[2]) = 0.0;
    Y(p[2]) = 0.0;
    Z(p[2]) = 0.0;
}

static inline void mw_addm(mwmatrix p, mwmatrix q, mwmatrix r)
{
    p[0] = mw_addv(q[0], r[0]);
    p[0] = mw_addv(q[1], r[1]);
    p[0] = mw_addv(q[2], r[2]);
}

static inline void mw_incaddm(mwmatrix p, mwmatrix q)
{
    mw_incaddv(p[0], q[0]);
    mw_incaddv(p[1], q[1]);
    mw_incaddv(p[2], q[2]);
}

static inline void mw_subm(mwmatrix p, mwmatrix q, mwmatrix r)
{
    p[0] = mw_subv(q[0], r[0]);
    p[0] = mw_subv(q[1], r[1]);
    p[0] = mw_subv(q[2], r[2]);
}

static inline void mw_incsubm(mwmatrix p, mwmatrix q)
{
    mw_incsubv(p[0], q[0]);
    mw_incsubv(p[1], q[1]);
    mw_incsubv(p[2], q[2]);
}

/* MULtiply Matrix by Scalar */
static inline void mw_mulms(mwmatrix p, mwmatrix q, real s)
{
    p[0] = mw_mulvs(q[0], s);
    p[1] = mw_mulvs(q[1], s);
    p[2] = mw_mulvs(q[2], s);
}

static inline void mw_incmulms(mwmatrix p, real s)
{
    mw_incmulvs(p[0], s);
    mw_incmulvs(p[1], s);
    mw_incmulvs(p[2], s);
}


#endif /* _MILKYWAY_VECTORS_H_ */