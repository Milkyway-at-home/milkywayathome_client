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

typedef struct MW_ALIGN_TYPE_V(4*sizeof(real))
{
    real x;
    real y;
    real z;
    real w;
} mwvector;

#define L(v) ((v)->x)
#define B(v) ((v)->y)
#define R(v) ((v)->z)

#define X(v) ((v)->x)
#define Y(v) ((v)->y)
#define Z(v) ((v)->z)
#define W(v) ((v)->w)

/* Cylinderical Macros */
#define CR(v) ((v)->x)
#define CT(v) ((v)->y)
#define CZ(v) ((v)->z)

#define SET_VECTOR(v, x, y, z) { X(v) = (x); Y(v) = (y); Z(v) = (z); W(v) = ZERO_REAL; }

#define mw_vec(x, y, z) { x, y, z, ZERO_REAL }


#define NDIM 3

#define ZERO_VECTOR { ZERO_REAL, ZERO_REAL, ZERO_REAL, ZERO_REAL }
typedef mwvector mwmatrix[NDIM];
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

#define MWVECTOR_TYPE "Vector"


CONST_F
static inline mwvector mw_addv(mwvector* a, mwvector* b)
{
    mwvector v;
    v.x = mw_add(&a->x, &b->x);
    v.y = mw_add(&a->y, &b->y);
    v.z = mw_add(&a->z, &b->z);

    return v;
}

CONST_F
static inline mwvector mw_subv(mwvector* a, mwvector* b)
{
    mwvector v;
    v.x = mw_sub(&a->x, &b->x);
    v.y = mw_sub(&a->y, &b->y);
    v.z = mw_sub(&a->z, &b->z);

    return v;
}

CONST_F
static inline mwvector mw_mulv(mwvector* a, mwvector* b)
{
    mwvector v;
    v.x = mw_mul(&a->x, &b->x);
    v.y = mw_mul(&a->y, &b->y);
    v.z = mw_mul(&a->z, &b->z);

    return v;
}

CONST_F
static inline mwvector mw_divv(mwvector* a, mwvector* b)
{
    mwvector v;
    v.x = mw_div(&a->x, &b->x);
    v.y = mw_div(&a->y, &b->y);
    v.z = mw_div(&a->z, &b->z);

    return v;
}

CONST_F
static inline real mw_dotv(mwvector* a, mwvector* b)
{
    real result = ZERO_REAL;
    result = mw_mad(&a->x, &b->x, &result);
    result = mw_mad(&a->y, &b->y, &result);
    result = mw_mad(&a->z, &b->z, &result);

    return result;
}

CONST_F
static inline mwvector mw_crossv(mwvector* a, mwvector* b)
{
    mwvector tmp;
    real bzay = mw_mul(&b->z, &a->y);
    real byaz = mw_mul(&b->y, &a->z);
    real bxaz = mw_mul(&b->x, &a->z);
    real bzax = mw_mul(&b->z, &a->x);
    real byax = mw_mul(&b->y, &a->x);
    real bxay = mw_mul(&b->x, &a->y);

    tmp.x = mw_sub(&bzay, &byaz);
    tmp.y = mw_sub(&bxaz, &bzax);
    tmp.z = mw_sub(&byax, &bxay);
    return tmp;
}

CONST_F
static inline real mw_length(mwvector* a)
{
    real R = mw_hypot(&a->x, &a->y);
    return mw_hypot(&R, &a->z);
}

CONST_F
static inline real mw_sqrv(mwvector* a)
{
    return mw_dotv(a, a);
}

CONST_F
static inline real mw_absv(mwvector* a)
{
    real sqrvec = mw_sqrv(a);
    return mw_sqrt(&sqrvec);
}

CONST_F
static inline mwvector mw_mulvs(mwvector* a, real* s)
{
    mwvector v;
    v.x =  mw_mul(&X(a), s);
    v.y =  mw_mul(&Y(a), s);
    v.z =  mw_mul(&Z(a), s);

    return v;
}

CONST_F
static inline mwvector mw_divvs(mwvector* a, real* s)
{
    mwvector v;
    v.x =  mw_div(&X(a), s);
    v.y =  mw_div(&Y(a), s);
    v.z =  mw_div(&Z(a), s);

    return v;
}

CONST_F
static inline mwvector mw_negv(mwvector* a)
{
    mwvector v;
    v.x = mw_neg(&a->x);
    v.y = mw_neg(&a->y);
    v.z = mw_neg(&a->z);

    return v;
}

CONST_F
static inline mwvector mw_mulmv(const mwmatrix m, mwvector* a)
{
    mwvector tmp;
    mwvector mx = m[0];
    mwvector my = m[1];
    mwvector mz = m[2];

    tmp.x = mw_dotv(&mx, a);
    tmp.y = mw_dotv(&my, a);
    tmp.z = mw_dotv(&mz, a);
    return tmp;
}

CONST_F
static inline real mw_distv(mwvector* a, mwvector* b)
{
    mwvector diff = mw_subv(a,b);
    real dist = mw_length(&diff);

    return dist;
}

/* Angle between two vectors, in the range [0,pi] */
CONST_F
static inline real mw_vecangle(mwvector* a, mwvector* b)
{
    real anorm, bnorm, dot;

    anorm = mw_length(a);
    bnorm = mw_length(b);
    dot = mw_dotv(a, b);

    real AB = mw_mul(&anorm, &bnorm);
    real costheta = mw_div(&dot, &AB);

    return mw_acos(&costheta);
}

CONST_F
static inline void mw_incsubv(mwvector* v1, mwvector* v2)
{
    *v1 = mw_subv(v1, v2);
}

CONST_F
static inline void mw_incaddv(mwvector* v1, mwvector* v2)
{
    *v1 = mw_addv(v1, v2);
}

CONST_F
static inline void mw_incmulv(mwvector* v1, mwvector* v2)
{
    *v1 = mw_mulv(v1, v2);
}

CONST_F
static inline void mw_incdivv(mwvector* v1, mwvector* v2)
{
    *v1 = mw_divv(v1, v2);
}

CONST_F
static inline void mw_incnegv(mwvector* v1)
{
    *v1 = mw_negv(v1);
}

#define mw_zerov(v) { (v)->x = ZERO_REAL; (v)->y = ZERO_REAL; (v)->z = ZERO_REAL; }

/* v1 += s * v2 */
CONST_F
static inline void mw_incaddv_s(mwvector* v1, mwvector* v2, real* s)
{
    mwvector tempor = mw_mulvs(v2, s);
    mw_incaddv(v1, &tempor);
}

/* v1 -= s * v2 */
CONST_F
static inline void mw_incsubv_s(mwvector* v1, mwvector* v2, real* s)
{
    mwvector tempor = mw_mulvs(v2, s);
    mw_incsubv(v1, &tempor);
}

CONST_F
static inline void mw_incdivs(mwvector* v, real* s)
{
    v->x = mw_div(&v->x, s);
    v->y = mw_div(&v->y, s);
    v->z = mw_div(&v->z, s);
}

CONST_F
static inline void mw_incmulvs(mwvector* v, real* s)
{
    v->x = mw_mul(&v->x, s);
    v->y = mw_mul(&v->y, s);
    v->z = mw_mul(&v->z, s);
}

CONST_F
static inline void mw_normalize(mwvector* v)
{
    real len = mw_length(v);
    mw_incdivs(v, &len);
}

/* Outer product */

static inline void mw_outv(mwmatrix p, mwvector* v, mwvector* u)
{
    X(&p[0]) = mw_mul(&X(v), &X(u));
    Y(&p[0]) = mw_mul(&X(v), &Y(u));
    Z(&p[0]) = mw_mul(&X(v), &Z(u));

    X(&p[1]) = mw_mul(&Y(v), &X(u));
    Y(&p[1]) = mw_mul(&Y(v), &Y(u));
    Z(&p[1]) = mw_mul(&Y(v), &Z(u));

    X(&p[2]) = mw_mul(&Z(v), &X(u));
    Y(&p[2]) = mw_mul(&Z(v), &Y(u));
    Z(&p[2]) = mw_mul(&Z(v), &Z(u));
}

/* Outer product of a vector with itself */
static inline void mw_outsqrv(mwmatrix p, mwvector* v)
{
    X(&p[0]) = mw_mul(&X(v), &X(v));
    Y(&p[0]) = mw_mul(&X(v), &Y(v));
    Z(&p[0]) = mw_mul(&X(v), &Z(v));

    X(&p[1]) = Y(&p[0]);
    Y(&p[1]) = mw_mul(&Y(v), &Y(v));
    Z(&p[1]) = mw_mul(&Y(v), &Z(v));

    X(&p[2]) = Z(&p[0]);
    Y(&p[2]) = Z(&p[1]);
    Z(&p[2]) = mw_mul(&Z(v), &Z(v));
}

static inline void mw_set_diagonal_matrix(mwmatrix p, real* s)
{
    X(&p[0]) = *s;
    Y(&p[0]) = ZERO_REAL;
    Z(&p[0]) = ZERO_REAL;

    X(&p[1]) = ZERO_REAL;
    Y(&p[1]) = *s;
    Z(&p[1]) = ZERO_REAL;

    X(&p[2]) = ZERO_REAL;
    Y(&p[2]) = ZERO_REAL;
    Z(&p[2]) = *s;
}

static inline void mw_set_matrix_identity(mwmatrix p)
{
    X(&p[0]) = mw_real_const(1.0);
    Y(&p[0]) = ZERO_REAL;
    Z(&p[0]) = ZERO_REAL;

    X(&p[1]) = ZERO_REAL;
    Y(&p[1]) = mw_real_const(1.0);
    Z(&p[1]) = ZERO_REAL;

    X(&p[2]) = ZERO_REAL;
    Y(&p[2]) = ZERO_REAL;
    Z(&p[2]) = mw_real_const(1.0);
}

static inline void mw_set_matrix_zero(mwmatrix p)
{
    X(&p[0]) = ZERO_REAL;
    Y(&p[0]) = ZERO_REAL;
    Z(&p[0]) = ZERO_REAL;

    X(&p[1]) = ZERO_REAL;
    Y(&p[1]) = ZERO_REAL;
    Z(&p[1]) = ZERO_REAL;

    X(&p[2]) = ZERO_REAL;
    Y(&p[2]) = ZERO_REAL;
    Z(&p[2]) = ZERO_REAL;
}

static inline void mw_addm(mwmatrix p, mwmatrix q, mwmatrix r)
{
    p[0] = mw_addv(&q[0], &r[0]);
    p[1] = mw_addv(&q[1], &r[1]);
    p[2] = mw_addv(&q[2], &r[2]);
}

static inline void mw_incaddm(mwmatrix p, mwmatrix q)
{
    mw_incaddv(&p[0], &q[0]);
    mw_incaddv(&p[1], &q[1]);
    mw_incaddv(&p[2], &q[2]);
}

static inline void mw_subm(mwmatrix p, mwmatrix q, mwmatrix r)
{
    p[0] = mw_subv(&q[0], &r[0]);
    p[0] = mw_subv(&q[1], &r[1]);
    p[0] = mw_subv(&q[2], &r[2]);
}

static inline void mw_incsubm(mwmatrix p, mwmatrix q)
{
    mw_incsubv(&p[0], &q[0]);
    mw_incsubv(&p[1], &q[1]);
    mw_incsubv(&p[2], &q[2]);
}

/* MULtiply Matrix by Scalar */
static inline void mw_mulms(mwmatrix p, mwmatrix q, real* s)
{
    p[0] = mw_mulvs(&q[0], s);
    p[1] = mw_mulvs(&q[1], s);
    p[2] = mw_mulvs(&q[2], s);
}

static inline void mw_incmulms(mwmatrix p, real* s)
{
    mw_incmulvs(&p[0], s);
    mw_incmulvs(&p[1], s);
    mw_incmulvs(&p[2], s);
}


#endif /* _MILKYWAY_VECTORS_H_ */
