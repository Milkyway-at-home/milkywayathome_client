/* ************************************************************************** */
/* MILKYWAY_VECTORS.H: include file for vector/matrix operations. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
   Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
   Rensselaer Polytechnic Institute. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _MILKYWAY_VECTORS_CPU_H_
#define _MILKYWAY_VECTORS_CPU_H_

#include "real.h"

#define NDIM 3

#define VECTOR_SIZE 3
typedef real vector[NDIM];

typedef real matrix[NDIM][NDIM];
typedef real* vectorptr;

/*
typedef real _mwvector __attribute__((vector_size(sizeof(real) * 4)));

typedef union
{
    real __attribute__ ((aligned(sizeof(real) * 4))) s[4];
    struct { real x, y, z, w; };
    _mwvector v4;
} mwvector;
*/


#define VECTOR(x, y, z) { (x), (y), (z) }

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
        (v)[0] = 0.0;                           \
        (v)[1] = 0.0;                           \
        (v)[2] = 0.0;                           \
    }



#endif /* _MILKYWAY_VECTORS_CPU_H_ */

