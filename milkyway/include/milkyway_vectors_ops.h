/* ************************************************************************** */
/* MILKYWAY_VECTORS.H: include file for vector/matrix operations. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
   Rensselaer Polytechnic Institute. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _MILKYWAY_VECTORS_OPS_H_
#define _MILKYWAY_VECTORS_OPS_H_

#include "real.h"

/* Vector operations. */

/* UNIT Vector */
#define UNITV(v,j)                              \
    {                                           \
        X(v) = (0 == (j)) ? 1.0 : 0.0;          \
        Y(v) = (1 == (j)) ? 1.0 : 0.0;          \
        Z(v) = (2 == (j)) ? 1.0 : 0.0;          \
    }

/* SET Vector */
#define SETV(v,u)                               \
    {                                           \
        X(v) = X(u);                            \
        Y(v) = Y(u);                            \
        Z(v) = Z(u);                            \
    }

/* ADD Vector */
#define ADDV(v,u,w)                             \
    {                                           \
        X(v) = X(u) + X(w);                     \
        Y(v) = Y(u) + Y(w);                     \
        Z(v) = Z(u) + Z(w);                     \
    }

/* SUBtract Vector */
#define SUBV(v,u,w)                             \
    {                                           \
        X(v) = X(u) - X(w);                     \
        Y(v) = Y(u) - Y(w);                     \
        Z(v) = Z(u) - Z(w);                     \
    }

/* INCSUBtract Vector  from vector multiplied scalar,
   v -= s * w;
*/
#define INCSUBVMS(v,s,u)                        \
    {                                           \
        X(v) -= (s) * X(u);                     \
        Y(v) -= (s) * Y(u);                     \
        Z(v) -= (s) * Z(u);                     \
    }


/* MULtiply Vector by Scalar */
#define MULVS(v,u,s)                            \
    {                                           \
        X(v) = (s) * X(u);                      \
        Y(v) = (s) * Y(u);                      \
        Z(v) = (s) * Z(u);                      \
    }

/* Negate vector */
#define NEGV(v,u)                               \
    {                                           \
        X(v) = -X(u);                           \
        Y(v) = -Y(u);                           \
        Z(v) = -Z(u);                           \
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
        X(v) = X(u) / (s);                  \
        Y(v) = Y(u) / (s);                  \
        Z(v) = Z(u) / (s);                  \
    }

/* DOT Vector Product */
#define DOTVP(s,v,u)                                    \
    {                                                   \
        (s) = X(v) * X(u) + Y(v) * Y(u) + Z(v) * Z(u);  \
    }

/* DOT Vector Product with itself */
#define SQRV(s,v)                                   \
    {                                               \
        (s) = sqr(X(v)) + sqr(Y(v)) + sqr(Z(v));    \
    }

/* ABSolute value of a Vector */
#define ABSV(s,v)                               \
    {                                           \
        real _tmp;                              \
        _tmp = sqr(X(v));                       \
        _tmp += sqr(Y(v));                      \
        _tmp += sqr(Z(v));                      \
        (s) = rsqrt(_tmp);                      \
    }

/* DISTance between Vectors */
#define DISTV(s,u,v)                            \
    {                                           \
        real _tmp;                              \
        _tmp = sqr(X(u)-X(v));                  \
        _tmp += sqr(Y(u)-Y(v));                 \
        _tmp += sqr(Z(u)-Z(v));                 \
        (s) = rsqrt(_tmp);                      \
    }

/* CROSS Vector Product */
#define CROSSVP(v,u,w)                          \
    {                                           \
        X(v) = Y(u)*Z(w) - Z(u)*Y(w);           \
        Y(v) = Z(u)*X(w) - X(u)*Z(w);           \
        Z(v) = X(u)*Y(w) - Y(u)*X(w);           \
    }

/* INCrementally ADD Vector */
#define INCADDV(v,u)                            \
    {                                           \
        X(v) += X(u);                           \
        Y(v) += Y(u);                           \
        Z(v) += Z(u);                           \
    }

/* INCrementally SUBtract Vector */
#define INCSUBV(v,u)                            \
    {                                           \
        X(v) -= X(u);                           \
        Y(v) -= Y(u);                           \
        Z(v) -= Z(u);                           \
    }

/* INCrementally MULtiply Vector by Scalar */
#define INCMULVS(v,s)                           \
    {                                           \
        X(v) *= (s);                            \
        Y(v) *= (s);                            \
        Z(v) *= (s);                            \
    }

/* INCrementally DIVide Vector by Scalar */
#define INCDIVVS(v,s)                           \
    {                                           \
        X(v) /= (s);                            \
        Y(v) /= (s);                            \
        Z(v) /= (s);                            \
    }

/* Matrix operations. */

#define CLRM(p)         /* CLeaR Matrix */      \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = 0.0;              \
    }

/* SET Matrix to Identity */
#define SETMI(p)                                        \
    {                                                   \
        size_t _i, _j;                                  \
        for (_i = 0; _i < NDIM; ++_i)                   \
            for (_j = 0; _j < NDIM; ++_j)               \
                (p)[_i][_j] = (_i == _j) ? 1.0 : 0.0;   \
    }

/* SET Matrix */
#define SETM(p,q)                               \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = (q)[_i][_j];      \
    }

/* TRANspose Matrix */
#define TRANM(p,q)                              \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = (q)[_j][_i];      \
    }

/* ADD Matrix */
#define ADDM(p,q,r)                                         \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i)                       \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (p)[_i][_j] = (q)[_i][_j] + (r)[_i][_j];    \
    }

/* SUBtract Matrix */
#define SUBM(p,q,r)                                         \
    {                                                       \
        size_t _i, _j;                                      \
        for (_i = 0; _i < NDIM; ++_i)                       \
            for (_j = 0; _j < NDIM; ++_j)                   \
                (p)[_i][_j] = (q)[_i][_j] - (r)[_i][_j];    \
    }

/* Multiply Matrix */
#define MULM(p,q,r)                                             \
    {                                                           \
        size_t _i, _j, _k;                                      \
        for (_i = 0; _i < NDIM; ++_i)                           \
            for (_j = 0; _j < NDIM; ++_j) {                     \
                (p)[_i][_j] = 0.0;                              \
                for (_k = 0; _k < NDIM; _k++)                   \
                    (p)[_i][_j] += (q)[_i][_k] * (r)[_k][_j];   \
            }                                                   \
    }

/* MULtiply Matrix by Scalar */
#define MULMS(p,q,s)                                \
    {                                               \
        size_t _i, _j;                              \
        for (_i = 0; _i < NDIM; ++_i)               \
            for (_j = 0; _j < NDIM; ++_j)           \
                (p)[_i][_j] = (q)[_i][_j] * (s);    \
    }

/* DIVide Matrix by Scalar */
#define DIVMS(p,q,s)                                \
    {                                               \
        size_t _i, _j;                              \
        for (_i = 0; _i < NDIM; ++_i)               \
            for (_j = 0; _j < NDIM; ++_j)           \
                (p)[_i][_j] = (q)[_i][_j] / (s);    \
    }

/* MULtiply Matrix by Vector */
#define MULMV(v,p,u)                                \
    {                                               \
        size_t _i, _j;                              \
        for (_i = 0; _i < NDIM; ++_i) {             \
            (v)[_i] = 0.0;                          \
            for (_j = 0; _j < NDIM; ++_j)           \
                (v)[_i] += (p)[_i][_j] * (u)[_j];   \
        }                                           \
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
        X(v) = (s);                             \
        Y(v) = (s);                             \
        Z(v) = (s);                             \
    }

/* ADD Vector and Scalar */
#define ADDVS(v,u,s)                            \
    {                                           \
        X(v) = X(u) + (s);                      \
        Y(v) = Y(u) + (s);                      \
        Z(v) = Z(u) + (s);                      \
    }

/* SET Matrix to Scalar */
#define SETMS(p,s)                              \
    {                                           \
        size_t _i, _j;                          \
        for (_i = 0; _i < NDIM; ++_i)           \
            for (_j = 0; _j < NDIM; ++_j)       \
                (p)[_i][_j] = (s);              \
    }

#endif /* _MILKYWAY_VECTORS_OPS_H_ */

