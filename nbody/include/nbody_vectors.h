/* ************************************************************************** */
/* MILKYWAY_VECTORS.H: include file for vector/matrix operations. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
   Rensselaer Polytechnic Institute. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _NBODY_VECTORS_H_
#define _NBODY_VECTORS_H_

/* FIXME: This entire thing should be killed with fire and replaced
 * with mwvec stuff */

#undef X
#undef Y
#undef Z
#undef W
#undef L
#undef B
#undef R
#undef ZERO_VECTOR
#undef ZERO_MATRIX


#define NDIM 3

#define VECTOR_SIZE 3
typedef real vector[NDIM];

typedef real matrix[NDIM][NDIM];
typedef real* vectorptr;


#define VECTOR(x, y, z) { (x), (y), (z) }
#define SET_VECTOR(v, x, y, z) { X(v) = (x); Y(v) = (y); Z(v) = (z); }

/* Note: vectorptr is NOT the same as vector*.  By using real* as
   vectorptr, we can do nice things to avoid pointer aliasing and
   copying in various places. Use vector* only for mapping over an
   array of vectors.
 */

#define ZERO_VECTOR { 0.0, 0.0, 0.0 }
#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

#define L(x) ((x)[0])
#define B(x) ((x)[1])
#define R(x) ((x)[2])

#define X(x) ((x)[0])
#define Y(x) ((x)[1])
#define Z(x) ((x)[2])


/* Clear vector */
#define CLRV(v)                                 \
    {                                           \
        X(v) = 0.0;                             \
        Y(v) = 0.0;                             \
        Z(v) = 0.0;                             \
    }

/* DOT Vector Product */
#define DOTVP(s,v,u)                                    \
    {                                                   \
        (s) = X(v) * X(u) + Y(v) * Y(u) + Z(v) * Z(u);  \
    }


/* Vector operations. */

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
        X(v) = -X(v);                           \
        Y(v) = -Y(v);                           \
        Z(v) = -Z(v);                           \
    }

/* DIVide Vector by Scalar */
#define DIVVS(v,u,s)                            \
    {                                           \
        X(v) = X(u) / (s);                      \
        Y(v) = Y(u) / (s);                      \
        Z(v) = Z(u) / (s);                      \
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
        (s) = mw_sqrt(_tmp);                    \
    }

/* DISTance between Vectors */
#define DISTV(s,u,v)                            \
    {                                           \
        real _tmp;                              \
        _tmp = sqr(X(u)-X(v));                  \
        _tmp += sqr(Y(u)-Y(v));                 \
        _tmp += sqr(Z(u)-Z(v));                 \
        (s) = mw_sqrt(_tmp);                    \
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

/* CLeaR Matrix */
#define CLRM(p)                                 \
    {                                           \
        X((p)[0]) = 0.0;                        \
        Y((p)[0]) = 0.0;                        \
        Z((p)[0]) = 0.0;                        \
                                                \
        X((p)[1]) = 0.0;                        \
        Y((p)[1]) = 0.0;                        \
        Z((p)[1]) = 0.0;                        \
                                                \
        X((p)[2]) = 0.0;                        \
        Y((p)[2]) = 0.0;                        \
        Z((p)[2]) = 0.0;                        \
    }

/* SET Matrix to Identity */
#define SETMI(p)                                \
    {                                           \
        X((p)[0]) = 1.0;                        \
        Y((p)[0]) = 0.0;                        \
        Z((p)[0]) = 0.0;                        \
                                                \
        X((p)[1]) = 0.0;                        \
        Y((p)[1]) = 1.0;                        \
        Z((p)[1]) = 0.0;                        \
                                                \
        X((p)[2]) = 0.0;                        \
        Y((p)[2]) = 0.0;                        \
        Z((p)[2]) = 1.0;                        \
    }

/* ADD Matrix */
#define ADDM(p,q,r)                                 \
    {                                               \
        size_t _i;                                  \
        for (_i = 0; _i < NDIM; ++_i)               \
        {                                           \
            X((p)[_i]) = X((q)[_i]) + X((r)[_i]);   \
            Y((p)[_i]) = Y((q)[_i]) + Y((r)[_i]);   \
            Z((p)[_i]) = Z((q)[_i]) + Z((r)[_i]);   \
        }                                           \
    }


/* SUBtract Matrix */


#define SUBM(p,q,r)                                 \
    {                                               \
        size_t _i;                                  \
        for (_i = 0; _i < NDIM; ++_i)               \
        {                                           \
            X((p)[_i]) = X((q)[_i]) - X((r)[_i]);   \
            Y((p)[_i]) = Y((q)[_i]) - Y((r)[_i]);   \
            Z((p)[_i]) = Z((q)[_i]) - Z((r)[_i]);   \
        }                                           \
    }

/* MULtiply Matrix by Scalar */
#define MULMS(p,q,s)                            \
    {                                           \
        size_t _i;                              \
        for (_i = 0; _i < NDIM; ++_i)           \
        {                                       \
            X((p)[_i]) = s * X((q)[_i]);        \
            Y((p)[_i]) = s * Y((q)[_i]);        \
            Z((p)[_i]) = s * Z((q)[_i]);        \
        }                                       \
    }

#define OUTVP(p,v,u)                            \
    {                                           \
        X((p)[0]) = X(v) * X(u);                \
        Y((p)[0]) = X(v) * Y(u);                \
        Z((p)[0]) = X(v) * Z(u);                \
                                                \
        X((p)[1]) = Y(v) * X(u);                \
        Y((p)[1]) = Y(v) * Y(u);                \
        Z((p)[1]) = Y(v) * Z(u);                \
                                                \
        X((p)[2]) = Z(v) * X(u);                \
        Y((p)[2]) = Z(v) * Y(u);                \
        Z((p)[2]) = Z(v) * Z(u);                \
    }


/* MULtiply Matrix by Vector */
#define MULMV(v, p, u)                          \
    {                                           \
        DOTVP(X(v), (p)[0], (u));               \
        DOTVP(Y(v), (p)[1], (u));               \
        DOTVP(Z(v), (p)[2], (u));               \
   }

#endif /* _NBODY_VECTORS_H_ */

