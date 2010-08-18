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

#include "real.h"

#define NDIM 3
#define VECTOR_SIZE 4

#ifdef __OPENCL_VERSION__ /* In the kernel */

  #if DOUBLEPREC
    typedef double4 real4, *real4ptr;
  #else
    typedef float4 real4, *real4ptr;
  #endif /* DOUBLEPREC */

  #define ZERO_VECTOR (0.0, 0.0, 0.0, 0.0)
  #define VECTOR(x, y, z) ( (x), (y), (z), 0.0 )

  #define L(x) ((x).x)
  #define B(x) ((x).y)
  #define R(x) ((x).z)

  #define X(x) ((x).x)
  #define Y(x) ((x).y)
  #define Z(x) ((x).z)

#else  /* Host */

  #if DOUBLEPREC
    typedef cl_double4 real4, *real4ptr;
  #else
    typedef cl_float4 real4, *real4ptr;
  #endif /* DOUBLEPREC */

  #define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
  #define VECTOR(x, y, z) { (x), (y), (z), 0.0 }

  #define L(x) ((x)[0])
  #define B(x) ((x)[1])
  #define R(x) ((x)[2])

  #define X(x) ((x)[0])
  #define Y(x) ((x)[1])
  #define Z(x) ((x)[2])

#endif /* __OPENCL_VERSION__ */

typedef real4 vector;
typedef real* vectorptr;

typedef real4 matrix[NDIM];

#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

/* clear vector. The 4th component can be ignored for everything else.

 */
#define CLRV(v)                                 \
    {                                           \
        (v)[0] = 0.0;                           \
        (v)[1] = 0.0;                           \
        (v)[2] = 0.0;                           \
        (v)[2] = 0.0;                           \
    }

#endif /* _MILKYWAY_VECTORS_CL_H_ */

