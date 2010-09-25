/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#else
  #include <time.h>
  #include <sys/time.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "milkyway_config.h"
#include "milkyway_util.h"
#include "milkyway_math.h"

#ifdef __SSE__
  #include <xmmintrin.h>
#endif /* __SSE__ */


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

#ifndef _WIN32
void* mwMallocAligned(size_t size, size_t alignment)
{
    void* p;

    if (posix_memalign(&p, alignment, size))
    {
        perror(__func__);
        fail("Failed to allocate block of size %zu aligned to %zu\n", size, alignment);
    }

    if (!p)
        fail("%s: NULL\n", __func__);

    return p;
}
#else
void* mwMallocAligned(size_t size, size_t alignment)
{
    void* p;
    p = _aligned_malloc(size, alignment);

    if (!p)
        fail("%s: NULL: _aligned_malloc error = %d\n", FUNC_NAME, GetLastError());

    return p;
}
#endif

void* mwCallocAligned(size_t count, size_t size, size_t alignment)
{
    void* p;

    p = mwMallocAligned(count * size, alignment);
    memset(p, 0, size);

    return p;
}

void* reallocSafe(void* ptr, size_t size)
{
    void* mem = (void*) realloc(ptr, size);
    if (mem == NULL)
        fail("realloc failed: %lu bytes\n", (unsigned long) size);
    return mem;
}

char* mwReadFile(const char* filename)
{
    FILE* f;
    long fsize;
    size_t readSize;
    char* buf;

    f = mw_fopen(filename, "rb");
    if (!f)
    {
        warn("Failed to open file '%s' for reading\n", filename);
        return NULL;
    }

    fseek(f, 0, SEEK_END);  /* Find size of file */
    fsize = ftell(f);

    fseek(f, 0, SEEK_SET);

    buf = callocSafe(fsize + 1, sizeof(char));

    readSize = fread(buf, sizeof(char), fsize, f);

    if (readSize != fsize)
    {
        free(buf);
        warn("Failed to read file '%s': Expected to read %ld, but got %u\n",
             filename,
             fsize,
             (unsigned int) readSize);
        return NULL;
    }

    return buf;
}

#if BOINC_APPLICATION

FILE* mwOpenResolved(const char* filename, const char* mode)
{
    int ret;
    char resolvedPath[1024];

    ret = boinc_resolve_filename(filename, resolvedPath, sizeof(resolvedPath));
    if (ret)
    {
        warn("Error resolving file '%s': %d\n", filename, ret);
        return NULL;
    }

    return mw_fopen(resolvedPath, mode);
}

#else

FILE* mwOpenResolved(const char* filename, const char* mode)
{
    return mw_fopen(filename, mode);
}

#endif /* BOINC_APPLICATION */


#ifdef _WIN32

double mwGetTime()
{
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return (double)t.QuadPart/(double)f.QuadPart;
}

#else

double mwGetTime()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

#endif

/* Modified from boinc_rename, which doesn't use MoveFileEx on
 * windows, which is more atomic.*/
static inline int _mwRename(const char* oldf, const char* newf)
{
  #ifdef _WIN32
    if (MoveFileEx(oldf, newf, MOVEFILE_REPLACE_EXISTING))
        return 0;
    return GetLastError();
  #else
    return rename(oldf, newf);
  #endif
}


#if BOINC_APPLICATION

int mwRename(const char* oldf, const char* newf)
{
    int rc;
    unsigned int i;

    /* FIXME: BOINC has random timing for retries. Fix boinc rename on
     * windows, then we can just get rid of this. */
    rc = _mwRename(oldf, newf);
    if (rc)
    {
        for (i = 0; i < 5; ++i)
        {
            sleep(1);       /* sleep 1 second, avoid lockstep */
            rc = _mwRename(oldf, newf);
            if (!rc)
                break;
        }
    }

    return rc;
}

#else

int mwRename(const char* oldf, const char* newf)
{
    return _mwRename(oldf, newf);
}

#endif /* BOINC_APPLICATION*/

#if defined(__SSE__) && DISABLE_DENORMALS

int mwDisableDenormalsSSE()
{
    int oldMXCSR = _mm_getcsr();
    int newMXCSR = oldMXCSR | 0x8040;
    _mm_setcsr(newMXCSR);

    warn("Disabled denormals\n");
    return oldMXCSR;
}

#endif /* defined(__SSE__) && DISABLE_DENORMALS */

/* From the extra parameters, read them as doubles */
real* mwReadRestArgs(const char** rest, const unsigned int numParams, unsigned int* pCountOut)
{
    unsigned int i;
    real* parameters = NULL;
    unsigned int paramCount = 0;

    /* Read through all the server arguments, and make sure we can
     * read everything and have the right number before trying to
     * do anything with them */

    if (numParams == 0)
    {
        warn("numParams = 0 makes no sense\n");
        return NULL;
    }

    if (!rest)
    {
        warn("%s: got rest == NULL\n", FUNC_NAME);
        return NULL;
    }

    while (rest[++paramCount]);  /* Count number of parameters */

    /* Make sure the number of extra parameters matches the number
     * we were told to expect. */
    if (numParams != paramCount)
    {
        warn("Parameter count mismatch: Expected %u, got %u\n", numParams, paramCount);
        return NULL;
    }

    parameters = (real*) mallocSafe(sizeof(real) * numParams);

    errno = 0;
    for (i = 0; i < numParams; ++i)
    {
        parameters[i] = (real) strtod(rest[i], NULL);
        if (errno)
        {
            perror("Error parsing command line fit parameters");
            free(parameters);
            return NULL;
        }
    }

    if (pCountOut)
        *pCountOut = paramCount;

    return parameters;
}

void _mw_time_prefix(char* buf, size_t bufSize)
{
    time_t x = time(NULL);
    struct tm* tm = localtime(&x);
    if (!strftime(buf, bufSize - 1, "%H:%M:%S", tm))
        buf[0] = '\0';
}

