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

#include "milkyway_config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifndef _WIN32
  #include <unistd.h>
  #include <fcntl.h>
#else
  #define _CRT_SECURE_NO_WARNINGS
  #define WIN32_LEAN_AND_MEAN
  #define VC_EXTRALEAN
  #include <malloc.h>
  #include <windows.h>
#endif /* _WIN32 */

#include <popt.h>

#include "milkyway_extra.h"
#include "milkyway_math.h"
#include "milkyway_show.h"
#include "mw_boinc_util.h"
#include "dSFMT.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MWInvalidEnum (-1)
#define InvalidEnum MWInvalidEnum


#if MILKYWAY_OPENCL

typedef struct
{
    cl_device_type devType;
    cl_uint platform;
    cl_uint devNum;
    cl_bool nonResponsive;  /* If screen redraws aren't important. Either don't care or something like an outputless Tesla */
    cl_uint numChunk;
    cl_double responsivenessFactor;
    cl_double targetFrequency;
    cl_int pollingMode;
    cl_bool enableCheckpointing;

    cl_bool forceNoIntrinsics;
    cl_bool forceX87;
    cl_bool forceSSE2;
    cl_bool forceSSE3;
    cl_bool verbose;
} CLRequest;

#else

/* CAL or nothing */
typedef struct
{
    unsigned int devNum;
    int nonResponsive;
    double responsivenessFactor;
    double targetFrequency;
    int pollingMode;
    int enableCheckpointing;

    int forceNoIntrinsics;
    int forceX87;
    int forceSSE2;
    int forceSSE3;
    int verbose;
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
  #if defined(_MSC_VER) || defined(__MINGW64__)
    #define mwFreeA _aligned_free
  #elif defined(__MINGW32__)
    #define mwFreeA __mingw_aligned_free
  #endif
#endif /* _WIN32 */


#define warn(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__)

/* Have value of 1 so we can do a return warn("blah blah\n") */
#define warn1(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__), 1

/* Controlled, but lazy failure */
#define fail(msg, ...)                              \
    {                                               \
        fprintf(stderr, msg, ##__VA_ARGS__);        \
        mw_finish(EXIT_FAILURE);                    \
    }

/* Failure related to a known code limitation */
#define mw_panic(msg, ...)                                              \
    {                                                                   \
        fprintf(stderr, "PANIC: in function '%s' %s(%d): " msg,         \
                FUNC_NAME, __FILE__, __LINE__, ##__VA_ARGS__);          \
        mw_finish(EXIT_FAILURE);                                        \
    }

#define mw_unreachable()                                                               \
    {                                                                                  \
        fprintf(stderr, "PANIC: Unreachable point reached in function '%s' %s(%d)\n",  \
                FUNC_NAME, __FILE__, __LINE__);                                        \
        mw_finish(EXIT_FAILURE);                                                       \
    }



/* If one of these options is null, use the default. */
#define stringDefault(s, d) ((s) = (s) ? (s) : strdup((d)))

char* mwReadFile(const char* filename);
char* mwFreadFile(FILE* f, const char* filename);
int mwWriteFile(const char* filename, const char* str);

double mwGetTime();
double mwGetTimeMilli();

#ifdef _WIN32
  #define mwMilliSleep(x) Sleep((x))
#else
#define mwMilliSleep(x) usleep(1000 * (x))
#endif /* _WIN32 */

int mwSetTimerMinResolution();
int mwResetTimerResolution();

#ifndef _WIN32

/* Just use nice value */
typedef int MWPriority;

#define MW_PRIORITY_DEFAULT 0

#else

/* The numbers these correspond to aren't usable */
typedef enum
{
    MW_PRIORITY_IDLE         = 0,
    MW_PRIORITY_BELOW_NORMAL = 1,
    MW_PRIORITY_NORMAL       = 2,
    MW_PRIORITY_ABOVE_NORMAL = 3,
    MW_PRIORITY_HIGH         = 4
} MWPriority;

#define MW_PRIORITY_DEFAULT MW_PRIORITY_NORMAL

#endif /* _WIN32 */

/* Default priority in case of invalid priority */
#ifndef _WIN32
#else

#endif /* _WIN32 */


int mwSetProcessPriority(MWPriority priority);

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
  #define mw_win_perror(str) fprintf(stderr, str ": %ld\n", GetLastError());
#endif

/* mw_xrandom: generate floating-point random number */
#define mwXrandom(st, xl, xh) ((real) (xl) + (real) ((xh) - (xl)) * dsfmt_genrand_open_open((st)))
#define mwUnitRandom(st) mwXrandom(st, -1.0, 1.0)

mwvector mwRandomUnitVector(dsfmt_t* dsfmtState);
mwvector mwRandomVector(dsfmt_t* dsfmtState, real r);

mwvector mwRandomUnitPoint(dsfmt_t* dsfmtState);
mwvector mwRandomPoint(dsfmt_t* dsfmtState, real r);


size_t mwDivRoundup(size_t a, size_t b);
#define mwEven(x) ((x) % 2 == 0)
#define mwDivisible(x, n) ((x) % (n) == 0)


int mwCheckNormalPosNum(real n);
int mwCheckNormalPosNumEps(real n);

const char** mwFixArgv(int argc, const char** argv);

/* Loop through all arguments and report bad arguments */
int mwReadArguments(poptContext context);

/* Read array of strings into doubles. Returns NULL on failure. */
real* mwReadRestArgs(const char** rest,            /* String array as returned by poptGetArgs() */
                     const unsigned int numParams, /* Expected number of parameters */
                     unsigned int* paramCountOut); /* (Optional) return count of actual number parameters that could have been read */

const char** mwGetForwardedArguments(const char** args, unsigned int* nForwardedArgs);

#if defined(__SSE__) && DISABLE_DENORMALS
int mwDisableDenormalsSSE();
#endif /* defined(__SSE__) && DISABLE_DENORMALS */

void mwSetConsistentx87FPUPrecision();

#ifdef __cplusplus
}
#endif


#endif /* _MILKYWAY_UTIL_H_ */

