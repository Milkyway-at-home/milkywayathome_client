/* ************************************************************************** */
/* MILKYWAY_VECTORS.H: include file for vector/matrix operations. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
   Rensselaer Polytechnic Institute. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _MILKYWAY_VECTORS_H_
#define _MILKYWAY_VECTORS_H_

#include "real.h"

#define NDIM 3

#if NBODY_OPENCL || SEPARATION_OPENCL || defined(__OPENCL_VERSION__)
  #define VECTOR_SIZE 4

  #ifdef __OPENCL_VERSION__ /* In the kernel */
    #if DOUBLEPREC
      typedef double4 real4, *real4ptr;
    #else
      typedef float4 real4, *real4ptr;
    #endif /* DOUBLEPREC */
  #else
    #if DOUBLEPREC
      typedef cl_double4 real4, *real4ptr;
    #else
      typedef cl_float4 real4, *real4ptr;
    #endif /* DOUBLEPREC */
  #endif /* __OPENCL_VERSION__ */

  typedef real4 vector;
  typedef real* vectorptr;

  typedef real4 matrix[NDIM];
  #define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
  #define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }
#else
  #define VECTOR_SIZE 3
  typedef real vector[NDIM];

  typedef real matrix[NDIM][NDIM];
  typedef real* vectorptr;

/* Note: vectorptr is NOT the same as vector*.  By using real* as
   vectorptr, we can do nice things to avoid pointer aliasing and
   copying in various places. Use vector* only for mapping over an
   array of vectors.
 */

  #define ZERO_VECTOR { 0.0, 0.0, 0.0 }
  #define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }
#endif /* NBODY_OPENCL */

#define L(x) (((vectorptr) (x))[0])
#define B(x) (((vectorptr) (x))[1])
#define R(x) (((vectorptr) (x))[2])

#define X(x) (((vectorptr) (x))[0])
#define Y(x) (((vectorptr) (x))[1])
#define Z(x) (((vectorptr) (x))[2])


/* Vector operations. */

/* I think the extra 4th component is only an issue for clearing. The
 * rest of the time it's ignored. */
/* CLeaR Vector */
#if NBODY_OPENCL || SEPARATION_OPENCL || defined(__OPENCL_VERSION__)
  #define CLRV(v)                               \
    {                                           \
        (v)[0] = 0.0;                           \
        (v)[1] = 0.0;                           \
        (v)[2] = 0.0;                           \
        (v)[2] = 0.0;                           \
    }

#else
  #define CLRV(v)                               \
    {                                           \
        (v)[0] = 0.0;                           \
        (v)[1] = 0.0;                           \
        (v)[2] = 0.0;                           \
    }
#endif /* NBODY_OPENCL || defined(__OPENCL_VERSION__) */

/* UNIT Vector */
#define UNITV(v,j)                              \
    {                                           \
        (v)[0] = (0 == (j)) ? 1.0 : 0.0;        \
        (v)[1] = (1 == (j)) ? 1.0 : 0.0;        \
        (v)[2] = (2 == (j)) ? 1.0 : 0.0;        \
    }

/* SET Vector */
#define SETV(v,u)                               \
    {                                           \
        (v)[0] = (u)[0];                        \
        (v)[1] = (u)[1];                        \
        (v)[2] = (u)[2];                        \
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

/* SUBtract Vector */
#define SUBV(v,u,w)                             \
    {                                           \
        real* restrict _vp = (v);               \
        const real* restrict _up = (u);         \
        const real* restrict _wp = (w);         \
        *_vp++ = (*_up++) - (*_wp++);           \
        *_vp++ = (*_up++) - (*_wp++);           \
        *_vp   = (*_up) - (*_wp);               \
    }

/* MULtiply Vector by Scalar */
#define MULVS(v,u,s)                            \
    {                                           \
        real* restrict _vp = (v);               \
        const real* restrict _up = (u);         \
        *_vp++ = (*_up++) * (s);                \
        *_vp++ = (*_up++) * (s);                \
        *_vp   = (*_up) * (s);                  \
    }

/* Negate vector */
#define NEGV(v,u)                               \
    {                                           \
        real* restrict _vp = (v);               \
        const real* restrict _up = (u);         \
        *_vp++ = -(*_up++);                     \
        *_vp++ = -(*_up++);                     \
        *_vp   = -(*_up);                       \
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
        (v)[0] = (u)[0] / (s);                  \
        (v)[1] = (u)[1] / (s);                  \
        (v)[2] = (u)[2] / (s);                  \
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
        real _tmp;                              \
        _tmp = sqr((v)[0]);                     \
        _tmp += sqr((v)[1]);                    \
        _tmp += sqr((v)[2]);                    \
        (s) = rsqrt(_tmp);                      \
    }

/* DISTance between Vectors */
#define DISTV(s,u,v)                            \
    {                                           \
        real _tmp;                              \
        _tmp = sqr((u)[0]-(v)[0]);              \
        _tmp += sqr((u)[1]-(v)[1]);             \
        _tmp += sqr((u)[2]-(v)[2]);             \
        (s) = rsqrt(_tmp);                      \
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
        (v)[0] += (u)[0];                       \
        (v)[1] += (u)[1];                       \
        (v)[2] += (u)[2];                       \
    }

/* INCrementally SUBtract Vector */
#define INCSUBV(v,u)                            \
    {                                           \
        (v)[0] -= (u)[0];                       \
        (v)[1] -= (u)[1];                       \
        (v)[2] -= (u)[2];                       \
    }

/* INCrementally MULtiply Vector by Scalar */
#define INCMULVS(v,s)                           \
    {                                           \
        (v)[0] *= (s);                          \
        (v)[1] *= (s);                          \
        (v)[2] *= (s);                          \
    }

/* INCrementally DIVide Vector by Scalar */
#define INCDIVVS(v,s)                           \
    {                                           \
        (v)[0] /= (s);                          \
        (v)[1] /= (s);                          \
        (v)[2] /= (s);                          \
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

/* OUTer Vector Product */
#define OUTVP(p,v,u)                                \
    {                                               \
        size_t _i, _j;                              \
        for (_i = 0; _i < NDIM; ++_i)               \
            for (_j = 0; _j < NDIM; ++_j)           \
                (p)[_i][_j] = (v)[_i] * (u)[_j];    \
    }

/* TRACE of Matrix */
#define TRACEM(s,p)                             \
    {                                           \
        size_t _i;                              \
        (s) = 0.0;                              \
        for (_i = 0.0; _i < NDIM; ++_i)         \
            (s) += (p)[_i][_i];                 \
    }

/* Misc. impure operations. */

/* SET Vector to Scalar */
#define SETVS(v,s)                              \
    {                                           \
        (v)[0] = (s);                           \
        (v)[1] = (s);                           \
        (v)[2] = (s);                           \
    }

/* ADD Vector and Scalar */
#define ADDVS(v,u,s)                            \
    {                                           \
        (v)[0] = (u)[0] + (s);                  \
        (v)[1] = (u)[1] + (s);                  \
        (v)[2] = (u)[2] + (s);                  \
    }

/* SET Matrix to Scalar */
#define SETMS(p,s)                              \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = (s);              \
    }

#endif /* _MILKYWAY_VECTORS_H_ */

