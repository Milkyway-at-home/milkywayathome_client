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

#ifndef _MILKYWAY_UTIL_H_
#define _MILKYWAY_UTIL_H_

#ifndef _MSC_VER
  #define _GNU_SOURCE
#endif

#include "milkyway_config.h"
#include "milkyway_extra.h"
#include "milkyway_math.h"
#include "mw_boinc_util.h"
#include "dSFMT.h"

#ifndef _WIN32
  #include <unistd.h>
  #include <fcntl.h>
#else
  #define _CRT_SECURE_NO_WARNINGS
  #define WIN32_LEAN_AND_MEAN
  #define VC_EXTRALEAN
  #include <windows.h>
#endif /* _WIN32 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <popt.h>


#ifdef __cplusplus
extern "C" {
#endif


#if MILKYWAY_OPENCL

typedef struct
{
    cl_device_type devType;
    cl_uint platform;
    cl_uint devNum;
    cl_bool nonResponsive;  /* If screen redraws aren't important. Either don't care or something like an outputless Tesla */
    cl_uint numChunk;
} CLRequest;

#else

typedef struct
{
    int _useless;
} CLRequest;

#endif /* MILKYWAY_OPENCL */


#if MW_ENABLE_DEBUG
    /* convenient functions for printing debugging stuffs */
 #define mw_debug(msg, ...) fprintf(stderr, "%s():%d: ", FUNC_NAME, __LINE__); \
                             fprintf(stderr, msg, __VA_ARGS__);
 #define mw_debugmsg(msg) puts(msg)
#else
  #define mw_debug(msg, ...) ((void) 0)
  #define mw_debugmsg(msg, ...) ((void) 0)
#endif

/* Allocations with that abort everything on failure */
void* mwMalloc(size_t size);
void* mwCalloc(size_t count, size_t size);
void* mwRealloc(void* ptr, size_t size);

/* Safe allocations aligned to 16 */
void* mwMallocA(size_t size);
void* mwCallocA(size_t count, size_t size);

#ifndef _WIN32
  #define mwFreeA free
#else
  #define mwFreeA _aligned_free
#endif /* _WIN32 */

#define warn(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__)
#define fail(msg, ...) { fprintf(stderr, msg, ##__VA_ARGS__);  \
                         mw_finish(EXIT_FAILURE); }

/* If one of these options is null, use the default. */
#define stringDefault(s, d) ((s) = (s) ? (s) : strdup((d)))

char* mwReadFile(const char* filename);
char* mwFreadFile(FILE* f, const char* filename);

double mwGetTime();
double mwGetTimeMilli();

#ifndef _WIN32
long mwGetTimeMicro();
#endif /* _WIN32 */

#define mus_to_s(mus) ((double) (mus) / 1.0e6)

void _mw_time_prefix(char* buf, size_t bufSize);
#define mw_report(msg, ...)                             \
    {                                                   \
        char _buf[256];                                 \
        _mw_time_prefix(_buf, sizeof(_buf));            \
        fprintf(stderr, "%s: " msg, _buf, ##__VA_ARGS__);   \
    }

#ifndef _WIN32
  #define mw_win_perror perror
#else
  #define mw_win_perror(str) fprintf(stderr, str ": %d\n", GetLastError());
#endif

/* mw_xrandom: generate floating-point random number */
#define mwXrandom(st, xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * dsfmt_genrand_open_open((st)))
#define mwUnitRandom(st) mwXrandom(st, -1.0, 1.0)


const char** mwFixArgv(int argc, const char** argv);

/* Loop through all arguments and report bad arguments */
int mwReadArguments(poptContext context);

/* Read array of strings into doubles. Returns NULL on failure. */
real* mwReadRestArgs(const char** rest,            /* String array as returned by poptGetArgs() */
                     const unsigned int numParams, /* Expected number of parameters */
                     unsigned int* paramCountOut); /* (Optional) return count of actual number parameters that could have been read */

#if defined(__SSE__) && DISABLE_DENORMALS
int mwDisableDenormalsSSE();
#endif /* defined(__SSE__) && DISABLE_DENORMALS */

void mwSetConsistentx87FPUPrecision();

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_UTIL_H_ */

