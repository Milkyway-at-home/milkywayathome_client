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

#ifdef _WIN32
  #include <windows.h>
#else
  #include <sys/time.h>
  #include <sys/resource.h>
#endif


/* allocate: memory allocation with error checking. */
void* allocate(int nb)
{
    void* mem;

    mem = (void*) calloc(nb, 1);        /* calloc zeros memory */
    if (mem == NULL)
        fail("allocate: not enough memory (%d bytes)\n", nb);
    return mem;
}


/* Found on SO. No idea if the Windows atually works */
#ifdef _WIN32

double get_time()
{
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return double(t.QuadPart)/double(f.QuadPart);
}

#else

double get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

#endif

