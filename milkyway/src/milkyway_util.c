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

#ifndef _WIN32
  #include <sys/time.h>
#else
  #include <stdlib.h>
  #include <malloc.h>
  #include <windows.h>
#endif /* _WIN32 */

#include <time.h>
#include <errno.h>
#include <float.h>

#ifdef __SSE__
  #include <xmmintrin.h>
#endif /* __SSE__ */


#include "milkyway_util.h"
#include "mw_boinc_util.h"

void* mwCalloc(size_t count, size_t size)
{
    void* mem = (void*) calloc(count, size);
    if (mem == NULL)
        fail("calloc failed: "ZU" bytes\n", count * size);
    return mem;
}

void* mwMalloc(size_t size)
{
    void* mem = (void*) malloc(size);
    if (mem == NULL)
        fail("malloc failed: "ZU" bytes\n", size);
    return mem;
}


#ifndef __APPLE__

void* mwCallocA(size_t count, size_t size)
{
    void* p;
    size_t totalSize = count * size;

    p = mwMallocA(totalSize);
    memset(p, 0, totalSize);

    return p;
}

#else

void* mwCallocA(size_t count, size_t size)
{
    return mwCalloc(count, size);
}


#endif /* __APPLE__ */

#if defined(__APPLE__)

/* OS X already aligns everything to 16 bytes */
void* mwMallocA(size_t size)
{
    return mwMalloc(size);
}

#elif !defined(_WIN32)

void* mwMallocA(size_t size)
{
    void* p;

    if (posix_memalign(&p, 16, size))
    {
        perror(__func__);
        fail("Failed to allocate block of size %zu aligned to 16\n", size);
    }

    if (!p)
        fail("%s: NULL\n", __func__);

    return p;
}

#else

#if defined(__MINGW32__) && !defined(__MINGW64__)
  #define _aligned_malloc __mingw_aligned_malloc
#endif

void* mwMallocA(size_t size)
{
    void* p;

    p = _aligned_malloc(size, 16);
    if (!p)
        fail("%s: NULL: _aligned_malloc error = %ld\n", FUNC_NAME, GetLastError());

    return p;
}
#endif /* defined(__APPLE__) */

void* mwRealloc(void* ptr, size_t size)
{
    void* mem = (void*) realloc(ptr, size);
    if (mem == NULL)
        fail("realloc failed: "ZU" bytes\n", size);
    return mem;
}

static char* fcloseVerbose(FILE* f, const char* err)
{
    if (fclose(f))
        perror(err);

    return NULL;
}

char* mwFreadFile(FILE* f, const char* filename)
{
    long fsize;
    size_t readSize;
    char* buf;

    if (!f)
    {
        warn("Failed to open file '%s' for reading\n", filename);
        return NULL;
    }

     /* Find size of file */
    if (fseek(f, 0, SEEK_END) == -1)
        return fcloseVerbose(f, "Seeking file end");

    fsize = ftell(f);
    if (fsize == -1)
        return fcloseVerbose(f, "Getting file size");

    fseek(f, 0, SEEK_SET);

    buf = (char*) mwMalloc((fsize + 1) * sizeof(char));
    buf[fsize] = '\0';

    readSize = fread(buf, sizeof(char), fsize, f);
    if (readSize != (size_t) fsize)
    {
        warn("Failed to read file '%s': Expected to read %ld, but got "ZU"\n",
             filename, fsize, readSize);
        free(buf);
        buf = NULL;
    }

    fcloseVerbose(f, "Closing read file");

    return buf;
}

char* mwReadFile(const char* filename)
{
    return mwFreadFile(mw_fopen(filename, "rb"), filename);
}

int mwWriteFile(const char* filename, const char* str)
{
    FILE* f;
    int rc;

    f = mw_fopen(filename, "w");
    if (!f)
    {
        perror("Writing file");
        return 1;
    }

    rc = fputs(str, f);
    if (rc == EOF)
        warn("Error writing file '%s'\n", filename);

    fcloseVerbose(f, "Closing write file");
    return rc;
}


#ifdef _WIN32

double mwGetTime()
{
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return (double)t.QuadPart/(double)f.QuadPart;
}

double mwGetTimeMilli()
{
    return 1.0e3 * mwGetTime();
}

#else

/* Seconds */
double mwGetTime()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

/* Get time in microseconds */
long mwGetTimeMicro()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return 1000000 * t.tv_sec  + t.tv_usec;
}

/* Get time in milliseconds */
double mwGetTimeMilli()
{
    return (double) mwGetTimeMicro() / 1.0e3;
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


#if WINDOWS_USES_X87

void mwSetConsistentx87FPUPrecision()
{
#if defined(_WIN32) && defined(_MSC_VER)
    //unsigned int control_word_x87;
    //control_word_x87 = __control87_2(_PC_64
    /* Set x87 intermediate rounding precision to 64 bits */
    warn("Setting FPU precision\n");
    _controlfp(_MCW_PC, _PC_64);
#else
  #warning Setting FPU flags with MinGW broken
#endif
}

#else

void mwSetConsistentx87FPUPrecision() { }

#endif

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

    parameters = (real*) mwMalloc(sizeof(real) * numParams);

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


int mwReadArguments(poptContext context)
{
    int o;

    while ( ( o = poptGetNextOpt(context)) >= 0 );

    if (o < -1)
    {
        warn("Argument parsing error: %s: %s\n",
             poptBadOption(context, 0),
             poptStrerror(o));

        return 1;
    }

    return 0;
}


/* This is a really stupid workaround: The server sends double
 * arguments. Some of them happen to be negative numbers. Since
 * BOINC is stupid and appends extra arguments, e.g. --nthreads 4
 * or --device 0 for multithreaded or GPU applications to the
 * command line instead of prepending like it should do, you can't
 * use POPT_CONTEXT_POSIXMEHARDER. Popt tries to interpret
 * negative numbers as options when not using
 * POPT_CONTEXT_POSIXMEHARDER, which results in them being
 * interpreted as wrong options and then erroring.
 */

/* Horrible function to find the -p -np arguments, and take anything
 * after them and move them to the front */
const char** mwFixArgv(unsigned long argc, const char** argv)
{
    const char** argvCopy;
    const char** p;
    char* endP;

    unsigned long i, j;
    unsigned long np, remaining, appendedCount;
    unsigned long npCheck = 0;

    p = argv;
    while (*p && strncmp(*p, "-np", 3)) /* Find how many numbers we're expected to find */
        ++p, ++npCheck;

    if (!*p)  /* Probably not using the server arguments */
        return NULL;

    if (!*(++p)) /* Go to the actual value of np */
        return NULL;

    /* The next argument after np should be the number */
    np = strtoul(*p, &endP, 10);
    if (*endP != '\0')
    {
        perror("Reading np");
        return NULL;
    }

    /* -2: -1 from argv[0], -1 from -p arg, -1 from np value */
    remaining = (argc - 3) - npCheck;  /* All remaining arguments */
    if (np > remaining || np >= argc)  /* Careful of underflow */
        return NULL;

    ++p;   /* Move on to the p argument */
    if (*p && strncmp(*p, "-p", 2))
    {
        warn("Didn't find expected p argument\n");
        return NULL;
    }

    if (!*(++p))  /* Should be first actual number argument */
        return NULL;

    /* FIXME: Have no dependence on np. FIXME: BOINC really, really,
     * really shouldn't ever append arguments ever (should prepend)
     * and it should be fixed. */

    argvCopy = (const char**) mwCalloc(argc, sizeof(const char*));
    i = j = 0;  /* Index into copy argv, original argv */
    argvCopy[i++] = argv[j++];
    p += np;  /* Skip the arguments */

    appendedCount = remaining - np;  /* Extras added by BOINC */
    while (*p && i <= appendedCount)
        argvCopy[i++] = *p++;

    while (i < argc && j < argc)  /* Copy the rest of the arguments into an appropriate order */
        argvCopy[i++] = argv[j++];

    return argvCopy;
}

size_t mwDivRoundup(size_t a, size_t b)
{
    return (a % b != 0) ? a / b + 1 : a / b;
}

mwvector mwRandomVector(dsfmt_t* dsfmtState)
{
    /* pick from unit cube */
    mwvector vec;

    X(vec) = mwUnitRandom(dsfmtState);
    Y(vec) = mwUnitRandom(dsfmtState);
    Z(vec) = mwUnitRandom(dsfmtState);
    W(vec) = 0.0;

    return vec;
}

