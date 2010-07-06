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
    const real xp = X(r) + ctx->sunGCDist;

    L(lbR) = ratan2(Y(r), xp);
    B(lbR) = ratan2( Z(r), rsqrt( sqr(xp) + sqr(Y(r)) ) );
    R(lbR) = rsqrt(sqr(xp) + sqr(Y(r)) + sqr(Z(r)));

    if (L(lbR) < 0.0)
        L(lbR) += 2 * M_PI;

}

void cartesianToLbr(const NBodyCtx* ctx, vectorptr restrict lbR, const vectorptr restrict r)
{
    cartesianToLbr_rad(ctx, lbR, r);
    L(lbR) = r2d(L(lbR));
    B(lbR) = r2d(B(lbR));
}

inline static void _lbrToCartesian(vectorptr cart, const real l, const real b, const real r, const real sun)
{
    X(cart) = r * rcos(l) * rcos(b) - sun;
    Y(cart) = r * rsin(l) * rcos(b);
    Z(cart) = r * rsin(b);
}

void lbrToCartesian_rad(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr)
{
    _lbrToCartesian(cart, L(lbr), B(lbr), R(lbr), ctx->sunGCDist);
}

void lbrToCartesian(const NBodyCtx* ctx, vectorptr cart, const vectorptr lbr)
{
    _lbrToCartesian(cart, d2r(L(lbr)), d2r(B(lbr)), R(lbr), ctx->sunGCDist);
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
    return (double)t.QuadPart/(double)f.QuadPart;
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

