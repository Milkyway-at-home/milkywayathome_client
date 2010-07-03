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

void cartesianToLbr_rad(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    const real r0p = r[0] + ctx->sunGCDist;
    lbR[0] = ratan2(r[1], r0p);
    lbR[1] = ratan2(r[2], rsqrt(sqr(r0p) + sqr(r[1])));
    lbR[2] = rsqrt(sqr(r0p) + sqr(r[1]) + sqr(r[2]));

    if (lbR[0] < 0)
        lbR[0] += 2 * M_PI;
}

void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    cartesianToLbr_rad(ctx, lbR, r);
    L(lbR) = r2d(L(lbR));
    B(lbR) = r2d(B(lbR));
}

inline static void _lbrToCartesian(vectorptr cart, const real l, const real b, const real r, const real sun)
{
    cart[0] = r * rcos(l) * rcos(b) - sun;
    cart[1] = r * rsin(l) * rcos(b);
    cart[2] = r * rsin(b);
}

void lbrToCartesian_rad(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr)
{
    _lbrToCartesian(cart, lbr[1], lbr[2], lbr[0], ctx->sunGCDist);
}

void lbrToCartesian(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr)
{
    _lbrToCartesian(cart, d2r(lbr[1]), d2r(lbr[2]), lbr[0], ctx->sunGCDist);
}

void* callocSafe(size_t count, size_t size)
{
    void* mem = (void*) calloc(count, size);
    if (mem == NULL)
        fail("calloc failed: %zd bytes\n", count * size);
    return mem;
}

void* mallocSafe(size_t size)
{
    void* mem = (void*) malloc(size);
    if (mem == NULL)
        fail("malloc failed: %zd bytes\n", size);
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

