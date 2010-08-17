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
  #define WIN32_LEAN_AND_MEAN
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
        fail("calloc failed: %lu bytes\n", (unsigned long) count * size);
    return mem;
}

void* mallocSafe(size_t size)
{
    void* mem = (void*) malloc(size);
    if (mem == NULL)
        fail("malloc failed: %lu bytes\n", (unsigned long) size);
    return mem;
}

#if BOINC_APPLICATION

FILE* nbodyOpenResolved(const char* filename, const char* mode)
{
    int ret;
    char resolvedPath[1024];

    ret = boinc_resolve_filename(filename, resolvedPath, sizeof(resolvedPath));
    if (ret)
        fail("Error resolving file '%s': %d\n", filename, ret);

    return nbody_fopen(resolvedPath, mode);
}

#else

FILE* nbodyOpenResolved(const char* filename, const char* mode)
{
    return nbody_fopen(filename, mode);
}

#endif /* BOINC_APPLICATION */

char* nbodyReadFile(const char* filename)
{
    FILE* f;
    long fsize;
    size_t readSize;
    char* buf;

    f = nbody_fopen(filename, "rb");
    if (!f)
    {
        perror("nbodyReadFile nbody_fopen:");
        warn("Failed to open file '%s' for reading\n", filename);
        return NULL;
    }

    /* Find size of file */
    if (fseek(f, 0, SEEK_END) < 0)
    {
        perror("nbodyReadFile fseek end:");
        return NULL;
    }

    fsize = ftell(f);
    if (fsize == -1)
    {
        perror("nbodyReadFile ftell:");
        return NULL;
    }

    if (fseek(f, 0, SEEK_SET) < 0)
    {
        perror("nbodyReadFile fseek beginning:");
        return NULL;
    }

    buf = callocSafe(fsize + 1, sizeof(char));

    readSize = fread(buf, sizeof(char), fsize, f);
    if (readSize != (unsigned long) fsize)
    {
        perror("nbodyReadFile fread:");
        free(buf);
        fclose(f);
        warn("Failed to read file '%s': Expected to read %ld, but got %u\n",
             filename,
             fsize,
             (unsigned int) readSize);
        return NULL;
    }

    return buf;
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

