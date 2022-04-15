/*
 *  Copyright (c) 2010-2011 Matthew Arsenault
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MILKYWAY_UTIL_H_
#define _MILKYWAY_UTIL_H_

#include "milkyway_config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#if HAVE_UNISTD_H
  #include <unistd.h>
#endif

#if HAVE_FCNTL_H
  #include <fcntl.h>
#endif

#if HAVE_MALLOC_H
  #include <malloc.h>
#endif

#if HAVE_WINDOWS_H
  #include <windows.h>
#endif

#include <popt.h>

#include "milkyway_extra.h"
#include "milkyway_math.h"
#include "milkyway_show.h"
#include "milkyway_timing.h"
#include "milkyway_alloc.h"
#include "milkyway_boinc_util.h"
#include "milkyway_rename.h"
#include "milkyway_asprintf.h"
#include "dSFMT.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MWInvalidEnum (-1)
#define InvalidEnum MWInvalidEnum


typedef struct
{
    const char* preferredPlatformVendor;
    unsigned int platform;
    unsigned int devNum;

    int nonResponsive;    /* If screen redraws aren't important. Either don't care or something like an outputless Tesla */
    double targetFrequency;
    double gpuWaitFactor;
    int pollingMode;
    int enableCheckpointing;

    int forceNoOpenCL;
    int forceNoILKernel;
    int forceNoIntrinsics;
    int forceX87; /* FIXME: Not always x87. More like a weird "other" kind of thing */
    int forceSSE2;
    int forceSSE3;
    int forceSSE41;
    int forceAVX;
    int verbose;
    int enableProfiling;
} CLRequest;

#if MW_ENABLE_DEBUG
    /* convenient functions for printing debugging stuffs */
  #define mw_debug(msg, ...)                                        \
    do                                                              \
    {                                                               \
        fprintf(stderr, "%s():%d: ", FUNC_NAME, __LINE__);          \
        fprintf(stderr, msg, __VA_ARGS__);                          \
    } while (0)

  #define mw_debugmsg(msg) puts(msg)
#else
  #define mw_debug(msg, ...) ((void) 0)
  #define mw_debugmsg(msg, ...) ((void) 0)
#endif


/* BOINC likes things on stderr */
#define mw_printf(msg, ...) fprintf(stderr, msg, ##__VA_ARGS__)

void mwPerror(const char* fmt, ...);

#ifdef _WIN32
void mwPerrorW32(const char* fmt, ...);
#endif

/* Controlled, but lazy failure */
#define mw_fail(...)                                \
    do                                              \
    {                                               \
        fprintf(stderr, __VA_ARGS__);               \
        mw_finish(EXIT_FAILURE);                    \
    } while (0)

/* Failure related to a known code limitation */
#define mw_panic(msg, ...)                                              \
    do                                                                  \
    {                                                                   \
        fprintf(stderr, "PANIC: in function '%s' %s(%d): " msg,         \
                FUNC_NAME, __FILE__, __LINE__, ##__VA_ARGS__);          \
        mw_finish(EXIT_FAILURE);                                        \
    } while (0)

#define mw_unreachable()                                                               \
    do                                                                                 \
    {                                                                                  \
        fprintf(stderr, "PANIC: Unreachable point reached in function '%s' %s(%d)\n",  \
                FUNC_NAME, __FILE__, __LINE__);                                        \
        mw_finish(EXIT_FAILURE);                                                       \
    } while (0)



/* If one of these options is null, use the default. */
#define mwStringDefault(s, d) ((s) = (s) ? (s) : strdup((d)))

char* mwReadFile(const char* filename);
char* mwReadFileWithSize(const char* filename, size_t* sizeOut);

char* mwFreadFile(FILE* f, const char* filename);
char* mwFreadFileWithSize(FILE* f, const char* filename, size_t* sizeOut);

int mwWriteFile(const char* filename, const char* str);

size_t mwCountLinesInFile(FILE* f);


/* Polling modes for OpenCL */

/* After initial wait when using mwCLWaitForEvent() how long (in ms)
 * should it sleep before trying again when using MW_POLL_SLEEP_CL_WAIT_FOR_EVENTS
 */
#define MW_DEFAULT_POLLING_PERIOD 3

/* Use clWaitForEvents() for waiting unless we suspect the driver has the high CPU issue */
#define MW_POLL_WORKAROUND_CL_WAIT_FOR_EVENTS -2

/* Use clWaitForEvents() normally */
#define MW_POLL_CL_WAIT_FOR_EVENTS -1

/* Use clWaitForEvents() after using an initial sleep based on execution time estimate */
#define MW_POLL_SLEEP_CL_WAIT_FOR_EVENTS 0



/* The numbers these correspond to aren't usable */
typedef enum
{
    MW_PRIORITY_INVALID      = -1,
    MW_PRIORITY_IDLE         = 0,
    MW_PRIORITY_BELOW_NORMAL = 1,
    MW_PRIORITY_NORMAL       = 2,
    MW_PRIORITY_ABOVE_NORMAL = 3,
    MW_PRIORITY_HIGH         = 4
} MWPriority;

#define MW_PRIORITY_DEFAULT MW_PRIORITY_NORMAL


int mwSetProcessPriority(MWPriority priority);

#ifndef _WIN32
long mwGetTimeMicro(void);
#endif /* _WIN32 */

#define mus_to_s(mus) ((double) (mus) / 1.0e6)

void mwLocalTime(char* buf, size_t bufSize);
void mwLocalTimeFull(char* buf, size_t bufSize);

#define mw_report(msg, ...)                                 \
    do                                                      \
    {                                                       \
        char _buf[256];                                     \
        mwLocalTime(_buf, sizeof(_buf));                    \
        fprintf(stdout, "%s: " msg, _buf, ##__VA_ARGS__);   \
    } while (0)

/* mw_xrandom: generate floating-point random number */
#define mwXrandom(st, xl, xh) ((real_0) (xl) + ((real_0) (xh) - (real_0) (xl)) * dsfmt_genrand_open_open((st)))
#define mwUnitRandom(st) mwXrandom(st, -1.0, 1.0)

mwvector mwRandomUnitVector(dsfmt_t* dsfmtState);
mwvector mwRandomVector(dsfmt_t* dsfmtState, real* r);

mwvector mwRandomUnitPoint(dsfmt_t* dsfmtState);
mwvector mwRandomPoint(dsfmt_t* dsfmtState, real* r);

real constrainAngle(real* a);


#define mwEven(x) ((x) % 2 == 0)
#define mwDivisible(x, n) ((x) % (n) == 0)
#define mwIsPowerOfTwo(x) (((x) != 0) && (((x) & (~(x) + 1)) == (x)));

/* integer ceil(a / b) */
#define mwDivRoundup(a, b) (((a) % (b) != 0) ? (a) / (b) + 1 : (a) / (b))

/* Find next multiple of b that is >= n */
#define mwNextMultiple(b, n) (((n) % (b)) ? ((n) + ((b) - (n) % (b))) : (n))

#define mwMin(a, b) ((a) < (b) ? (a) : (b));
#define mwMax(a, b) ((a) < (b) ? (b) : (a));


int mwCheckNormalNum(real_0 n);
int mwCheckNormalPosNum(real_0 n);
int mwCheckNormalPosNumEps(real_0 n);

const char** mwFixArgv(int argc, const char* argv[]);

/* Loop through all arguments and report bad arguments */
int mwReadArguments(poptContext context);

/* Read array of strings into doubles. Returns NULL on failure. */
real_0* mwReadRestArgs(const char** rest, unsigned int n);

const char** mwGetForwardedArguments(const char** args, unsigned int* nForwardedArgs);

#if defined(__SSE__) && DISABLE_DENORMALS
int mwDisableDenormalsSSE(void);
#endif /* defined(__SSE__) && DISABLE_DENORMALS */

unsigned long long mwFixFPUPrecision(void);
void mwDisableErrorBoxes(void);

int mwProcessIsAlive(int pid);


#ifdef __GNUC__
  #define mw_likely(x)    __builtin_expect((x), 1)
  #define mw_unlikely(x)  __builtin_expect((x), 0)

  #if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)
    #define mw_assume_aligned(p, a) __builtin_assume_aligned(p, a)
  #else
    #define mw_assume_aligned(p, a) (p)
  #endif
#else
  #define mw_likely(x) (x)
  #define mw_unlikely(x) (x)
  #define mw_assume_aligned(p, a) (p)
#endif /* __GNUC__ */


#ifdef __cplusplus
}
#endif


#endif /* _MILKYWAY_UTIL_H_ */

