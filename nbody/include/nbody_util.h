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

real cputime(void);              /* return elapsed CPU time */
void* allocate(int);             /* allocate and zero memory */
void eprintf(char*, ...);        /* printf to error FILE* */
int compare (const void* a, const void* b);     /* comparison function used in chisq */
real chisq();                  /* likelihood calculator */

#endif /* _NBODY_UTIL_H_ */

