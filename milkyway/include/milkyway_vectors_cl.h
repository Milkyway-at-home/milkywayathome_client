/* ************************************************************************** */
/* MILKYWAY_VECTORS.H: include file for vector/matrix operations. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
   Rensselaer Polytechnic Institute. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _MILKYWAY_VECTORS_CL_H_
#define _MILKYWAY_VECTORS_CL_H_

#include "milkyway_cl.h"
#include "real.h"

#define NDIM 3
#define VECTOR_SIZE 4

#ifdef __OPENCL_VERSION__ /* In the kernel */

  #if DOUBLEPREC
    typedef double4 real4;
    typedef double2 real2;
  #else
    typedef float4 real4;
    typedef float2 real2;
  #endif /* DOUBLEPREC */

  #define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
  #define VECTOR(x, y, z) { (x), (y), (z), 0.0 }

/* DOT Vector Product */

#if !BROKEN_CL_MATH
  #define DOTVP(s,v,u)                                  \
    {                                                   \
        (s) = dot(v, u);                                \
    }
#else
  #warning "Using broken dot product"
  #define DOTVP(s,v,u)                                  \
    {                                                   \
        (s) = sqr(X(v)) + sqr(Y(v)) + sqr(Z(v));        \
    }
#endif /* !BROKEN_CL_MATH */

  #define L(v) ((v).x)
  #define B(v) ((v).y)
  #define R(v) ((v).z)

  #define X(v) ((v).x)
  #define Y(v) ((v).y)
  #define Z(v) ((v).z)
  #define W(v) ((v).w)


#else  /* Host */
  #if DOUBLEPREC
    typedef cl_double4 real4;
    typedef cl_double2 real2;
  #else
    typedef cl_float4 real4;
    typedef cl_float2 real2;
  #endif /* DOUBLEPREC */

  #define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
  #define VECTOR(x, y, z) { (x), (y), (z), 0.0 }

/* DOT Vector Product */
#define DOTVP(s,v,u)                                    \
    {                                                   \
        (s) = X(v) * X(u) + Y(v) * Y(u) + Z(v) * Z(u);  \
    }

  #ifdef __APPLE__
    /* The host side implementation of cl_double4 seems on OS X to be
     * just an array of cl_double. This seems to not be true on ATI's,
     * using a more complicated union of different structs, which is
     * actually sort of nicer */
    #define L(v) ((v)[0])
    #define B(v) ((v)[1])
    #define R(v) ((v)[2])

    #define X(v) ((v)[0])
    #define Y(v) ((v)[1])
    #define Z(v) ((v)[2])
    #define W(v) ((v)[3])

    /* clear vector. The 4th component can be ignored for everything else. */
    #define CLRV(v)                             \
    {                                           \
        (v)[0] = 0.0;                           \
        (v)[1] = 0.0;                           \
        (v)[2] = 0.0;                           \
        (v)[2] = 0.0;                           \
    }

  #else
    #define L(v) ((v).x)
    #define B(v) ((v).y)
    #define R(v) ((v).z)

    #define X(v) ((v).x)
    #define Y(v) ((v).y)
    #define Z(v) ((v).z)
    #define W(v) ((v).w)

    #define CLRV(v)                            \
    {                                          \
        (v).x = 0.0;                           \
        (v).y = 0.0;                           \
        (v).z = 0.0;                           \
        (v).w = 0.0;                           \
    }
  #endif /* __APPLE__ */

#endif /* __OPENCL_VERSION__ */

typedef real4 vector;
typedef real* vectorptr;

typedef real4 matrix[NDIM];

#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }


#endif /* _MILKYWAY_VECTORS_CL_H_ */

