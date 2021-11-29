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

#define SET_VECTOR(v, x, y, z) { X(v) = (x); Y(v) = (y); Z(v) = (z); }

#define mw_vec(x, y, z) { x, y, z, ZERO_REAL }


#define NDIM 3

#define ZERO_VECTOR { ZERO_REAL, ZERO_REAL, ZERO_REAL, ZERO_REAL }
typedef mwvector mwmatrix[NDIM];
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

#define MWVECTOR_TYPE "Vector"


CONST_F
static inline mwvector mw_addv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = mw_add(a.x, b.x);
    v.y = mw_add(a.y, b.y);
    v.z = mw_add(a.z, b.z);
    return v;
}

CONST_F
static inline mwvector mw_subv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = mw_sub(a.x, b.x);
    v.y = mw_sub(a.y, b.y);
    v.z = mw_sub(a.z, b.z);

    return v;
}

CONST_F
static inline mwvector mw_mulv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = mw_mul(a.x, b.x);
    v.y = mw_mul(a.y, b.y);
    v.z = mw_mul(a.z, b.z);

    return v;
}

CONST_F
static inline mwvector mw_divv(mwvector a, mwvector b)
{
    mwvector v;
    v.x = mw_div(a.x, b.x);
    v.y = mw_div(a.y, b.y);
    v.z = mw_div(a.z, b.z);

    return v;
}

CONST_F
static inline real mw_dotv(mwvector a, mwvector b)
{
    return mw_mad(a.z, b.z, mw_mad(a.y, b.y, mw_mul(a.x, b.x)));
}

CONST_F
static inline mwvector mw_crossv(mwvector a, mwvector b)
{
    mwvector tmp;
    tmp.x = mw_sub(mw_mul(b.z, a.y), mw_mul(b.y, a.z));
    tmp.y = mw_sub(mw_mul(b.x, a.z), mw_mul(b.z, a.x));
    tmp.z = mw_sub(mw_mul(b.y, a.x), mw_mul(b.x, a.y));
    return tmp;
}

CONST_F
static inline real mw_length(mwvector a)
{
    return mw_sqrt(mw_add(mw_add(sqr(a.x), sqr(a.y)), sqr(a.z)));
}

CONST_F
static inline real mw_sqrv(mwvector a)
{
    return mw_mad(a.z, a.z, mw_mad(a.y, a.y, sqr(a.x)));
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
    v.x =  mw_mul(s, a.x);
    v.y =  mw_mul(s, a.y);
    v.z =  mw_mul(s, a.z);
    return v;
}

CONST_F
static inline mwvector mw_divvs(mwvector a, real s)
{
    mwvector v;
    v.x = mw_div(a.x, s);
    v.y = mw_div(a.y, s);
    v.z = mw_div(a.z, s);
    return v;
}

CONST_F
static inline mwvector mw_negv(mwvector a)
{
    mwvector v;
    v.x = mw_mul(mw_real_const(-1.0),a.x);
    v.y = mw_mul(mw_real_const(-1.0),a.y);
    v.z = mw_mul(mw_real_const(-1.0),a.z);
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
    return mw_sqrt(mw_add(mw_add(sqr(mw_sub(u.x, v.x)), sqr(mw_sub(u.y, v.y))), sqr(mw_sub(u.z, v.z))));
}

/* Angle between two vectors, in the range [0,pi] */
CONST_F
static inline real mw_vecangle(mwvector a, mwvector b)
{
    real anorm, bnorm, dot;

    anorm = mw_length(a);
    bnorm = mw_length(b);
    dot = mw_dotv(a, b);

    return mw_acos(mw_div(dot, mw_mul(anorm, bnorm)));
}


#define mw_incsubv(v1, v2) { (v1).x = mw_sub((v1).x, (v2).x); (v1).y = mw_sub((v1).y, (v2).y); (v1).z = mw_sub((v1).z, (v2).z); }
#define mw_incaddv(v1, v2) { (v1).x = mw_add((v1).x, (v2).x); (v1).y = mw_add((v1).y, (v2).y); (v1).z = mw_add((v1).z, (v2).z); }
#define mw_incdivv(v1, v2) { (v1).x = mw_div((v1).x, (v2).x); (v1).y = mw_div((v1).y, (v2).y); (v1).z = mw_div((v1).z, (v2).z); }
#define mw_incmulv(v1, v2) { (v1).x = mw_mul((v1).x, (v2).x); (v1).y = mw_mul((v1).y, (v2).y); (v1).z = mw_mul((v1).z, (v2).z); }
#define mw_incnegv(v1) { (v1).x = mw_mul(mw_real_const(-1.0),(v1).x); (v1).y = mw_mul(mw_real_const(-1.0),(v1).y); (v1).z = mw_mul(mw_real_const(-1.0),(v1).z); }
#define mw_zerov(v) { (v).x = ZERO_REAL; (v).y = ZERO_REAL; (v).z = ZERO_REAL; }

/* v1 += s * v2 */
#define mw_incaddv_s(v1, v2, s) { mw_incaddv((v1),mw_mulvs((v2),(s))); }

/* v1 -= s * v2 */
#define mw_incsubv_s(v1, v2, s) { mw_incsubv((v1),mw_mulvs((v2),(s))); }

#define mw_incdivs(v, s) { (v).x = mw_div((v).x,(s)); (v).y = mw_div((v).y,(s)); (v).z = mw_div((v).z,(s)); }
#define mw_incmulvs(v, s) { (v).x = mw_mul((v).x,(s)); (v).y = mw_mul((v).y,(s)); (v).z = mw_mul((v).z,(s)); }

#define mw_normalize(v) { real len = mw_length(v); mw_incdivs(v, len); }

/* Outer product */

static inline void mw_outv(mwmatrix p, mwvector v, mwvector u)
{
    X(p[0]) = mw_mul(X(v), X(u));
    Y(p[0]) = mw_mul(X(v), Y(u));
    Z(p[0]) = mw_mul(X(v), Z(u));

    X(p[1]) = mw_mul(Y(v), X(u));
    Y(p[1]) = mw_mul(Y(v), Y(u));
    Z(p[1]) = mw_mul(Y(v), Z(u));

    X(p[2]) = mw_mul(Z(v), X(u));
    Y(p[2]) = mw_mul(Z(v), Y(u));
    Z(p[2]) = mw_mul(Z(v), Z(u));
}

/* Outer product of a vector with itself */
static inline void mw_outsqrv(mwmatrix p, mwvector v)
{
    X(p[0]) = mw_mul(X(v), X(v));
    Y(p[0]) = mw_mul(X(v), Y(v));
    Z(p[0]) = mw_mul(X(v), Z(v));

    X(p[1]) = Y(p[0]);
    Y(p[1]) = mw_mul(Y(v), Y(v));
    Z(p[1]) = mw_mul(Y(v), Z(v));

    X(p[2]) = Z(p[0]);
    Y(p[2]) = Z(p[1]);
    Z(p[2]) = mw_mul(Z(v), Z(v));
}

static inline void mw_set_diagonal_matrix(mwmatrix p, real s)
{
    X(p[0]) = s;
    Y(p[0]) = ZERO_REAL;
    Z(p[0]) = ZERO_REAL;

    X(p[1]) = ZERO_REAL;
    Y(p[1]) = s;
    Z(p[1]) = ZERO_REAL;

    X(p[2]) = ZERO_REAL;
    Y(p[2]) = ZERO_REAL;
    Z(p[2]) = s;
}

static inline void mw_set_matrix_identity(mwmatrix p)
{
    X(p[0]) = mw_real_const(1.0);
    Y(p[0]) = ZERO_REAL;
    Z(p[0]) = ZERO_REAL;

    X(p[1]) = ZERO_REAL;
    Y(p[1]) = mw_real_const(1.0);
    Z(p[1]) = ZERO_REAL;

    X(p[2]) = ZERO_REAL;
    Y(p[2]) = ZERO_REAL;
    Z(p[2]) = mw_real_const(1.0);
}

static inline void mw_set_matrix_zero(mwmatrix p)
{
    X(p[0]) = ZERO_REAL;
    Y(p[0]) = ZERO_REAL;
    Z(p[0]) = ZERO_REAL;

    X(p[1]) = ZERO_REAL;
    Y(p[1]) = ZERO_REAL;
    Z(p[1]) = ZERO_REAL;

    X(p[2]) = ZERO_REAL;
    Y(p[2]) = ZERO_REAL;
    Z(p[2]) = ZERO_REAL;
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
