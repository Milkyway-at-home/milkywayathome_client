/* ************************************************************************** */
/* UTIL: various useful routines and functions. */
/* */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI. */
/* It's free because it's yours. */
/* ************************************************************************** */

#include "real.h"
#include "nbody.h"
#include "nbody_util.h"
#include "vectmath.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

#ifndef HZ
#  include <time.h>
#  define HZ CLK_TCK
#endif

/* xrandom: generate floating-point random number. */
real xrandom(real xl, real xh)
{
    return xl + (xh - xl) * drand48();
}

/* cputime: compute CPU time in minutes. */
real cputime()
{
    struct tms buffer;

    if (times(&buffer) == -1)
        fail("times() call failed\n");
    return buffer.tms_utime / (60.0 * HZ);
}

/* allocate: memory allocation with error checking. */

void* allocate(int nb)
{
    void* mem;

    mem = (void*) calloc(nb, 1);        /* calloc zeros memory */
    if (mem == NULL)
        fail("allocate: not enough memory (%d bytes)\n", nb);
    return mem;
}


