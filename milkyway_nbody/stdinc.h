/* ************************************************************************** */
/* STDINC.H: standard include file for C programs. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _STDINC_H_
#define _STDINC_H_

typedef short int bool;

#ifndef TRUE
#  define TRUE  1
#  define FALSE 0
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* PROC, IPROC: pointers to procedures and integer-valued functions. */

typedef void (*proc)();
typedef int (*iproc)();

/*  ABS: returns the absolute value of its argument
 *  MAX: returns the argument with the highest value
 *  MIN: returns the argument with the lowest value
 */

#define   ABS(x)       (((x) < 0) ? -(x) : (x))
#define   MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ? (x) : (y))


/* degrees to radians */
#define d2r(x) (x * (M_PI / 180.0))

/* radians to degrees */
#define r2d(x) (x * (180.0 / M_PI))

/* simple math macros */
#define cube(x) (x * x * x)
#define sqr(x)  (x * x)
#define inv(x)  (1.0 / x)

/* other useful nonstandard constants */

/* (4 * pi) / 3 Should be enough for 128-bit long double should we choose to use
 * that */
#define PI_4_3 (4.1887902047863909846168578443727)

#endif /* _STDINC_H_ */

