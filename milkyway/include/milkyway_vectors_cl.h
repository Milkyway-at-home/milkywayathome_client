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

#ifndef __OPENCL_VERSION__
  #include <OpenCL/cl.h>
  #include <OpenCL/cl_ext.h>
#endif /* __OPENCL_VERSION__ */

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
  #define VECTOR(x, y, z) ((vector) ( (x), (y), (z), 0.0 ))

  #define L(v) ((v).x)
  #define B(v) ((v).y)
  #define R(v) ((v).z)

  #define X(v) ((v).x)
  #define Y(v) ((v).y)
  #define Z(v) ((v).z)

#else  /* Host */
  #if DOUBLEPREC
    typedef cl_double4 real4, *real4ptr;
  #else
    typedef cl_float4 real4, *real4ptr;
  #endif /* DOUBLEPREC */

  #define ZERO_VECTOR { 0.0, 0.0, 0.0, 0.0 }
  #define VECTOR(x, y, z) { (x), (y), (z), 0.0 }

  #define L(v) ((v)[0])
  #define B(v) ((v)[1])
  #define R(v) ((v)[2])

  #define X(v) ((v)[0])
  #define Y(v) ((v)[1])
  #define Z(v) ((v)[2])

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

