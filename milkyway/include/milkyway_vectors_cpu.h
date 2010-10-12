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



/* Note: vectorptr is NOT the same as vector*.  By using real* as
   vectorptr, we can do nice things to avoid pointer aliasing and
   copying in various places. Use vector* only for mapping over an
   array of vectors.
 */

#define ZERO_MATRIX { ZERO_VECTOR, ZERO_VECTOR, ZERO_VECTOR }

/* Clear vector */
#define CLRV(v)                                 \
    {                                           \
        (v)[0] = 0.0;                           \
        (v)[1] = 0.0;                           \
        (v)[2] = 0.0;                           \
    }




#endif /* _MILKYWAY_VECTORS_CPU_H_ */

