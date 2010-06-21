/* ************************************************************************** */
/* VECTMATH.H: include file for vector/matrix operations. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _VECTMATH_H_
#define _VECTMATH_H_

#include "real.h"

#define NDIM 3

typedef real vector[NDIM], matrix[NDIM][NDIM];

typedef real* vectorptr;

#define ZERO_VECTOR { 0.0, 0.0, 0.0 }
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

/* Vector operations. */

#define CLRV(v)         /* CLeaR Vector */      \
    {                                           \
        size_t _i;                              \
        for (_i = 0; _i < NDIM; ++_i)           \
            (v)[_i] = 0.0;                      \
    }

#define UNITV(v,j)      /* UNIT Vector */       \
    {                                           \
        size_t _i;                              \
        for (_i = 0; _i < NDIM; ++_i)           \
            (v)[_i] = (_i == (j)) ? 1.0 : 0.0;  \
    }

#define SETV(v,u)       /* SET Vector */        \
    {                                           \
        size_t _i;                              \
        for (_i = 0; _i < NDIM; ++_i)           \
            (v)[_i] = (u)[_i];                  \
    }

/* Warning: Don't add a vector to itself, and don't increment with
 * this */
#define ADDV(v,u,w)     /* ADD Vector */            \
    {                                               \
        real* restrict _vp = (v);                   \
        const real* restrict _up = (u);             \
        const real* restrict _wp = (w);             \
        *_vp++ = (*_up++) + (*_wp++);               \
        *_vp++ = (*_up++) + (*_wp++);               \
        *_vp   = (*_up) + (*_wp);                   \
    }

#define SUBV(v,u,w)     /* SUBtract Vector */                   \
    {                                                           \
        real* restrict _vp = (v);                               \
        const real* restrict _up = (u);                         \
        const real* restrict _wp = (w);                         \
        *_vp++ = (*_up++) - (*_wp++);                           \
        *_vp++ = (*_up++) - (*_wp++);                           \
        *_vp   = (*_up) - (*_wp);                               \
    }

#define MULVS(v,u,s)        /* MULtiply Vector by Scalar */     \
    {                                                           \
        real* restrict _vp = (v);                               \
        const real* restrict _up = (u);                         \
        *_vp++ = (*_up++) * (s);                                \
        *_vp++ = (*_up++) * (s);                                \
        *_vp   = (*_up) * (s);                                  \
    }

/* Negate vector */
#define NEGV(v,u)                                               \
    {                                                           \
        real* restrict _vp = (v);                               \
        const real* restrict _up = (u);                         \
        *_vp++ = -(*_up++);                                     \
        *_vp++ = -(*_up++);                                     \
        *_vp   = -(*_up);                                       \
    }

#define INCNEGV(v)                              \
    {                                           \
        v[0] = -v[0];                           \
        v[1] = -v[1];                           \
        v[2] = -v[2];                           \
    }

/* DIVide Vector by Scalar */
#define DIVVS(v,u,s)                            \
    {                                           \
        size_t _i;                              \
        for (_i = 0; _i < NDIM; ++_i)           \
            (v)[_i] = (u)[_i] / (s);            \
    }

/* DOT Vector Product */
/* Warning: Don't dot a vector with itself. Use SQRV */
#define DOTVP(s,v,u)                            \
    {                                           \
        real* restrict _vp = (v);               \
        const real* restrict _up = (u);         \
        (s)  = (*_vp++) * (*_up++);             \
        (s) += (*_vp++) * (*_up++);             \
        (s) += (*_vp)   * (*_up);               \
    }

/* DOT Vector Product with itself */
#define SQRV(s,v)                               \
    {                                           \
        const real* restrict _vp = (v);         \
        (s)  = (*_vp) * (*_vp);                 \
        ++_vp;                                  \
        (s) += (*_vp) * (*_vp);                 \
        ++_vp;                                  \
        (s) += (*_vp) * (*_vp);                 \
    }

/* ABSolute value of a Vector */
#define ABSV(s,v)                               \
    {                                           \
        size_t _i;                              \
        real _tmp = 0.0;                        \
        for (_i = 0; _i < NDIM; ++_i)           \
            _tmp += (v)[_i] * (v)[_i];          \
        (s) = rsqrt(_tmp);                      \
    }

/* DISTance between Vectors */
#define DISTV(s,u,v)                                        \
    {                                                       \
        size_t _i;                                          \
        real _tmp = 0.0;                                    \
        for (_i = 0; _i < NDIM; ++_i)                       \
            _tmp += ((u)[_i]-(v)[_i]) * ((u)[_i]-(v)[_i]);  \
        (s) = rsqrt(_tmp);                                  \
    }

/* CROSS Vector Product */
#define CROSSVP(v,u,w)                          \
    {                                           \
        (v)[0] = (u)[1]*(w)[2] - (u)[2]*(w)[1]; \
        (v)[1] = (u)[2]*(w)[0] - (u)[0]*(w)[2]; \
        (v)[2] = (u)[0]*(w)[1] - (u)[1]*(w)[0]; \
    }

/* INCrementally ADD Vector */
#define INCADDV(v,u)                            \
    {                                           \
        size_t _i;                              \
        for (_i = 0; _i < NDIM; ++_i)           \
            (v)[_i] += (u)[_i];                 \
    }

#define INCSUBV(v,u)             /* INCrementally SUBtract Vector */    \
    {                                                                   \
        size_t _i;                                                      \
        for (_i = 0; _i < NDIM; ++_i)                                   \
            (v)[_i] -= (u)[_i];                                         \
    }

#define INCMULVS(v,s)   /* INCrementally MULtiply Vector by Scalar */   \
    {                                                                   \
        size_t _i;                                                      \
        for (_i = 0; _i < NDIM; ++_i)                                   \
            (v)[_i] *= (s);                                             \
    }

#define INCDIVVS(v,s)   /* INCrementally DIVide Vector by Scalar */ \
    {                                                               \
        size_t _i;                                                  \
        for (_i = 0; _i < NDIM; ++_i)                               \
            (v)[_i] /= (s);                                         \
    }

/* Matrix operations. */

#define CLRM(p)         /* CLeaR Matrix */      \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = 0.0;              \
    }

#define SETMI(p)        /* SET Matrix to Identity */    \
    {                                                   \
        size_t _i, _j;                                  \
        for (_i = 0; _i < NDIM; ++_i)                   \
            for (_j = 0; _j < NDIM; ++_j)               \
                (p)[_i][_j] = (_i == _j) ? 1.0 : 0.0;   \
    }

#define SETM(p,q)       /* SET Matrix */        \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = (q)[_i][_j];      \
    }

#define TRANM(p,q)      /* TRANspose Matrix */  \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = (q)[_j][_i];      \
    }

#define ADDM(p,q,r)     /* ADD Matrix */                    \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i)                       \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (p)[_i][_j] = (q)[_i][_j] + (r)[_i][_j];    \
    }

#define SUBM(p,q,r)     /* SUBtract Matrix */               \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i)                       \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (p)[_i][_j] = (q)[_i][_j] - (r)[_i][_j];    \
    }

#define MULM(p,q,r)     /* Multiply Matrix */                   \
    {                                                           \
        size_t _i, _j, _k;                                      \
        for (_i = 0; _i < NDIM; ++_i)                           \
            for (_j = 0; _j < NDIM; ++_j) {                     \
                (p)[_i][_j] = 0.0;                              \
                for (_k = 0; _k < NDIM; _k++)                   \
                    (p)[_i][_j] += (q)[_i][_k] * (r)[_k][_j];   \
            }                                                   \
    }

#define MULMS(p,q,s)        /* MULtiply Matrix by Scalar */ \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i)                       \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (p)[_i][_j] = (q)[_i][_j] * (s);            \
    }

#define DIVMS(p,q,s)        /* DIVide Matrix by Scalar */   \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i)                       \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (p)[_i][_j] = (q)[_i][_j] / (s);            \
    }

#define MULMV(v,p,u)        /* MULtiply Matrix by Vector */ \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i) {                     \
            (v)[_i] = 0.0;                                  \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (v)[_i] += (p)[_i][_j] * (u)[_j];           \
        }                                                   \
    }

#define OUTVP(p,v,u)        /* OUTer Vector Product */  \
    {                                                   \
        size_t _i, _j;                                  \
        for (_i = 0; _i < NDIM; ++_i)                   \
            for (_j = 0; _j < NDIM; ++_j)               \
                (p)[_i][_j] = (v)[_i] * (u)[_j];        \
    }

#define TRACEM(s,p)     /* TRACE of Matrix */   \
    {                                           \
        size_t _i;                              \
        (s) = 0.0;                              \
        for (_i = 0.0; _i < NDIM; ++_i)         \
            (s) += (p)[_i][_i];                 \
    }

/*  * Misc. impure operations.
 */

#define SETVS(v,s)      /* SET Vector to Scalar */  \
    {                                               \
        size_t _i;                                  \
        for (_i = 0; _i < NDIM; ++_i)               \
            (v)[_i] = (s);                          \
    }

#define ADDVS(v,u,s)        /* ADD Vector and Scalar */     \
    {                                                       \
        size_t _i;                                          \
        for (_i = 0; _i < NDIM; ++_i)                       \
            (v)[_i] = (u)[_i] + (s);                        \
    }

#define SETMS(p,s)      /* SET Matrix to Scalar */  \
    {                                               \
        size_t _i, _j;                              \
        for (_i = 0; _i < NDIM; ++_i)               \
            for (_j = 0; _j < NDIM; ++_j)           \
                (p)[_i][_j] = (s);                  \
    }

#endif /* _VECTMATH_H_ */

