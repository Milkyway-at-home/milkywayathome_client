/* ************************************************************************** */
/* DEFS.H: include file for hierarchical force calculation routines.  The */
/* definitions in this file are needed for load.c and grav.c; this file */
/* does not provide definitions for other parts of the N-body code. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _NBODY_UTIL_H_
#define _NBODY_UTIL_H_

#include "stdinc.h"
#include "real.h"
#include "vectmath.h"

/* Utility routines used in load.c and grav.c.  These are defined in
 * util.c, which must be compiled with the same choice of precision.
 */

/* xrandom: generate floating-point random number. */
#define xrandom(xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * drand48())

void* callocSafe(size_t count, size_t size);
void* mallocSafe(size_t size);

__attribute__ ((visibility("default"))) double get_time();

#endif /* _NBODY_UTIL_H_ */

