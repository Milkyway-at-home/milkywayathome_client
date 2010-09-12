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

#include "milkyway_util.h"

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

char* mwReadFile(const char* filename)
{
    FILE* f;
    long fsize;
    size_t readSize;
    char* buf;

    f = mw_fopen(filename, "r");
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
        fail("Error resolving file '%s': %d\n", filename, ret);

    return mw_fopen(resolvedPath, mode);
}

#else

FILE* mwOpenResolved(const char* filename, const char* mode)
{
    return mw_fopen(filename, mode);
}

#endif /* BOINC_APPLICATION */


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



