/* ************************************************************************** */
/* STDINC.H: standard include file for C programs. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#ifndef _STDINC_H_
#define _STDINC_H_

#ifndef TRUE
typedef short int bool;
#  define FALSE 0
#  define TRUE  1
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

#include /* _STDINC_H_ */

