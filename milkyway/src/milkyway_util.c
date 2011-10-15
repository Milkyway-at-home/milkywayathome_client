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

#ifndef _WIN32
  #include <sys/time.h>
  #include <sys/resource.h>
#else
  #include <stdlib.h>
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


#if HAVE_MACH_ABSOLUTE_TIME
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#if MW_IS_X86
  #if HAVE_FPU_CONTROL_H
    #include <fpu_control.h>
  #endif
  #ifndef _FPU_SETCW
    #define _FPU_SETCW(cw) __asm__ ("fldcw %0" : : "m" (*&cw))
  #endif
  #ifndef _FPU_GETCW
    #define _FPU_GETCW(cw) __asm__ ("fnstcw %0" : "=m" (*&cw))
  #endif

  #ifndef _FPU_EXTENDED
    #define _FPU_EXTENDED 0x300
  #endif
  #ifndef _FPU_DOUBLE
    #define _FPU_DOUBLE 0x200
  #endif
  #if defined(__FreeBSD__) && !defined(_FPU_DEFAULT)
    #include <machine/fpu.h>
    #define _FPU_DEFAULT __INITIAL_FPUCW__
  #endif
#endif /* MW_IS_X86 */


static char* fcloseVerbose(FILE* f, const char* name, const char* err)
{
    if (fclose(f))
    {
        mwPerror("Error closing file '%s' while %s", name, err);
    }

    return NULL;
}

char* mwFreadFileWithSize(FILE* f, const char* filename, size_t* sizeOut)
{
    long fsize;
    size_t readSize;
    char* buf;

    if (!f)
    {
        return NULL;
    }

     /* Find size of file */
    if (fseek(f, 0, SEEK_END) == -1)
    {
        return fcloseVerbose(f, filename, "seeking file end");
    }

    fsize = ftell(f);
    if (fsize == -1)
    {
        return fcloseVerbose(f, filename, "getting file size");
    }

    fseek(f, 0, SEEK_SET);

    buf = (char*) mwMalloc((fsize + 1) * sizeof(char));
    buf[fsize] = '\0';

    readSize = fread(buf, sizeof(char), fsize, f);
    if (readSize != (size_t) fsize)
    {
        mwPerror("Failed to read file '%s': Expected to read %ld, but got "ZU,
                 filename, fsize, readSize);
        free(buf);
        buf = NULL;
    }

    fcloseVerbose(f, filename, "closing read file");

    if (sizeOut)
        *sizeOut = readSize;

    return buf;
}

char* mwFreadFile(FILE* f, const char* filename)
{
    return mwFreadFileWithSize(f, filename, NULL);
}

char* mwReadFileWithSize(const char* filename, size_t* sizeOut)
{
    return mwFreadFileWithSize(mw_fopen(filename, "rb"), filename, sizeOut);
}

char* mwReadFile(const char* filename)
{
    return mwFreadFileWithSize(mw_fopen(filename, "rb"), filename, NULL);
}


int mwWriteFile(const char* filename, const char* str)
{
    FILE* f;
    int rc;

    if (!str || !filename)
    {
        return 1;
    }

    f = mw_fopen(filename, "wb");
    if (!f)
    {
        mwPerror("Writing file '%s'", filename);
        return 1;
    }

    rc = fputs(str, f);
    if (rc == EOF)
    {
        mwPerror("Error writing file '%s'", filename);
    }

    fcloseVerbose(f, filename, "Closing write file");
    return rc;
}

size_t mwCountLinesInFile(FILE* f)
{
    int c;
    size_t lineCount = 0;

    while ((c = fgetc(f)) != EOF)
    {
        if (c == '\n')
        {
            ++lineCount;
        }
    }

    if (!feof(f))
    {
        mwPerror("Error counting file lines");
        return 0;
    }

    if (fseek(f, 0L, SEEK_SET) < 0)
    {
        mwPerror("Error seeking file for counting");
        return 0;
    }

    return lineCount;
}

#ifdef _WIN32

double mwGetTime(void)
{
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return (double)t.QuadPart/(double)f.QuadPart;
}

double mwGetTimeMilli(void)
{
    return 1.0e3 * mwGetTime();
}

static UINT mwGetMinTimerResolution(void)
{
    TIMECAPS tc;

    if (timeGetDevCaps(&tc, sizeof(tc)) != TIMERR_NOERROR)
    {
        mw_printf("Failed to get timer resolution\n");
        return 0;
    }

	/* target 1-millisecond target resolution */
    return MIN(MAX(tc.wPeriodMin, 1), tc.wPeriodMax);
}

int mwSetTimerMinResolution(void)
{
    if (timeBeginPeriod(mwGetMinTimerResolution()) != TIMERR_NOERROR)
    {
        mw_printf("Failed to set timer resolution\n");
        return 1;
    }

    return 0;
}

int mwResetTimerResolution(void)
{
    if (timeEndPeriod(mwGetMinTimerResolution()) != TIMERR_NOERROR)
    {
        mw_printf("Failed to end timer resolution\n");
        return 1;
    }

    return 0;
}

#else

/* Seconds */
double mwGetTime(void)
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    /* Prevent weird breakage when building with -fsingle-precision-constant */
    return t.tv_sec + t.tv_usec * (double) 1.0e-6;
}

/* Get time in microseconds */
long mwGetTimeMicro(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return 1000000 * t.tv_sec  + t.tv_usec;
}

/* Get time in milliseconds */
double mwGetTimeMilli(void)
{
    return (double) mwGetTimeMicro() / 1.0e3;
}

int mwSetTimerMinResolution(void)
{
	return 0;
}

int mwResetTimerResolution(void)
{
    return 0;
}

#endif

#if defined(__SSE__) && DISABLE_DENORMALS

int mwDisableDenormalsSSE(void)
{
    int oldMXCSR = _mm_getcsr();
    int newMXCSR = oldMXCSR | 0x8040;
    _mm_setcsr(newMXCSR);

    mw_printf("Disabled denormals\n");
    return oldMXCSR;
}

#endif /* defined(__SSE__) && DISABLE_DENORMALS */


/* From the extra parameters, read them as doubles */
real* mwReadRestArgs(const char** rest, unsigned int n)
{
    unsigned int i;
    real* parameters = NULL;

    if (!rest)
        return NULL;

    parameters = (real*) mwMalloc(n * sizeof(real));

    errno = 0;
    for (i = 0; i < n; ++i)
    {
        parameters[i] = (real) strtod(rest[i], NULL);
        if (errno)
        {
            mwPerror("Error parsing command line fit parameters at '%s'", rest[i]);
            free(parameters);
            return NULL;
        }
    }

    return parameters;
}

/* Take subset of argv not used by the application to forward on to the Lua scripts.
   Need to free return value.
 */
const char** mwGetForwardedArguments(const char** args, unsigned int* nForwardedArgs)
{
    unsigned int i, argCount = 0;
    const char** forwardedArgs;

    if (!args)
    {
        *nForwardedArgs = 0;
        return NULL;
    }

    while (args[++argCount]);  /* Count number of parameters */

    forwardedArgs = (const char**) mwMalloc(sizeof(const char*) * argCount);

    for (i = 0; i < argCount; ++i)
        forwardedArgs[i] = args[i];

    *nForwardedArgs = argCount;

    return forwardedArgs;
}

void mwLocalTime(char* buf, size_t bufSize)
{
    time_t x = time(NULL);
    struct tm* tm = localtime(&x);
    if (!strftime(buf, bufSize - 1, "%H:%M:%S", tm))
        buf[0] = '\0';
}

void mwLocalTimeFull(char* buf, size_t bufSize)
{
    time_t x = time(NULL);
    struct tm* tm = localtime(&x);
    if (!strftime(buf, bufSize - 1, "%c", tm))
        buf[0] = '\0';
}

/* Returns < 0 on error, otherwise bitwise or'd val flags from arguments encountered.
   We want to be able to some behaviour without the flag, but simply using some default value for the flag is undesirable.
*/
int mwReadArguments(poptContext context)
{
    int o;
    int rc = 0;

    while ((o = poptGetNextOpt(context)) >= 0)
    {
         rc |= o;
    }

    if (o < -1)
    {
        mw_printf("Argument parsing error: %s: %s\n",
                  poptBadOption(context, 0),
                  poptStrerror(o));
        return -1;
    }

    return rc;
}


/* This is a really stupid workaround: The server sends double
 * arguments. Some of them happen to be negative numbers. Since BOINC
 * is stupid and appends extra arguments, e.g. --nthreads 4 or
 * --device 0 for multithreaded or GPU applications to the command
 * line instead of prepending like it should do (also user arguments
 * in app_info.xml), you can't use POPT_CONTEXT_POSIXMEHARDER. Popt
 * tries to interpret negative numbers as options when not using
 * POPT_CONTEXT_POSIXMEHARDER, which results in them being interpreted
 * as wrong options and then erroring.
 */
const char** mwFixArgv(int argc, const char* argv[])
{
    const char** argvCopy = NULL;
    const char** p = NULL;
    char* endP = NULL;

    int i, j;
    int np, remaining, appendedCount;
    int npCheck = 0;

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
        mwPerror("Reading np");
        return NULL;
    }

    /* -2: -1 from argv[0], -1 from -p arg, -1 from np value */
    remaining = (argc - 3) - npCheck;  /* All remaining arguments */
    if (np > remaining || np >= argc)  /* Careful of underflow */
    {
        mw_printf("Warning: Number of parameters remaining can't match expected: -np = %d\n", np);
        return NULL;
    }

    ++p;   /* Move on to the p argument */
    if (*p && strncmp(*p, "-p", 2))
    {
        mw_printf("Didn't find expected p argument\n");
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

/* Pick from unit cube. i.e. not normalized */
mwvector mwRandomUnitPoint(dsfmt_t* dsfmtState)
{
    mwvector vec;

    X(vec) = mwUnitRandom(dsfmtState);
    Y(vec) = mwUnitRandom(dsfmtState);
    Z(vec) = mwUnitRandom(dsfmtState);
    W(vec) = 0.0;

    return vec;
}

mwvector mwRandomPoint(dsfmt_t* dsfmtState, real s)
{
    mwvector v = mwRandomUnitPoint(dsfmtState);
    mw_incmulvs(v, s);
    return v;
}

/* Random direction */
mwvector mwRandomUnitVector(dsfmt_t* dsfmtState)
{
    mwvector v = mwRandomUnitPoint(dsfmtState);
    mw_normalize(v);
    return v;
}

mwvector mwRandomVector(dsfmt_t* dsfmtState, real r)
{
    mwvector v = mwRandomUnitVector(dsfmtState);
    mw_incmulvs(v, r);
    return v;
}


/* Check for a timesteps etc. which will actually finish. */
int mwCheckNormalPosNumEps(real n)
{
    return !isfinite(n) || n <= 0.0 || n <= REAL_EPSILON;
}

/* Check for positive, real numbers that can be any size */
int mwCheckNormalPosNum(real n)
{
    return !isfinite(n) || n <= 0.0;
}

/* Disable the stupid "Windows is checking for a solution to the
 * problem" boxes from showing up if this crashes */
void mwDisableErrorBoxes(void)
{
  #ifdef _WIN32

    DWORD mode, setMode;

    setMode = SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX;
    mode = SetErrorMode(setMode);
    SetErrorMode(mode | setMode);
  #endif /* _WIN32 */
}


#ifdef _WIN32

static DWORD mwPriorityToPriorityClass(MWPriority x)
{
    switch (x)
    {
        case MW_PRIORITY_IDLE:
            return IDLE_PRIORITY_CLASS;
        case MW_PRIORITY_BELOW_NORMAL:
            return BELOW_NORMAL_PRIORITY_CLASS;
        case MW_PRIORITY_NORMAL:
            return NORMAL_PRIORITY_CLASS;
        case MW_PRIORITY_ABOVE_NORMAL:
            return ABOVE_NORMAL_PRIORITY_CLASS;
        case MW_PRIORITY_HIGH:
            return HIGH_PRIORITY_CLASS;
        default:
            mw_printf("Invalid priority: %d. Using default.\n", (int) x);
            return mwPriorityToPriorityClass(MW_PRIORITY_DEFAULT);
    }
}

int mwSetProcessPriority(MWPriority priority)
{
    HANDLE handle;
    DWORD pClass;

    if (!DuplicateHandle(GetCurrentProcess(),
                         GetCurrentProcess(),
                         GetCurrentProcess(),
                         &handle,
                         0,
                         FALSE,
                         DUPLICATE_SAME_ACCESS))
    {
        mwPerrorW32("Failed to get process handle");
        return 1;
    }

    pClass = mwPriorityToPriorityClass(priority);
    if (!SetPriorityClass(handle, pclass))
    {
        mwPerrorW32("Failed to set process priority class to %ld", pClass);
        return 1;
    }

    return 0;
}
#else

int mwSetProcessPriority(MWPriority priority)
{
    if (setpriority(PRIO_PROCESS, getpid(), priority))
    {
        mwPerror("Setting process priority to %d", priority);
        return 1;
    }

    return 0;
}

#endif /* _WIN32 */


/*  From crlibm */
unsigned long long mwFixFPUPrecision(void)
{
#if MW_IS_X86

  #if defined(_MSC_VER) || defined(__MINGW32__)
  unsigned int oldcw, cw;

  /* CHECKME */
  oldcw = _controlfp(0, 0);

  #ifdef __MINGW32__
  _controlfp(_PC_53, _MCW_PC);
  #else
  _controlfp(_PC_53, MCW_PC);
  #endif

  return (unsigned long long) oldcw;

  #elif defined(__SUNPRO_C)  /* Sun Studio  */
  unsigned short oldcw, cw;

  __asm__ ("movw    $639, -22(%ebp)");
  __asm__ ("fldcw -22(%ebp)");

  return (unsigned long long) oldcw;
  #elif !defined(__APPLE__) /* GCC, clang */
  /* x87 FPU never used on OS X */
  unsigned short oldcw, cw;

    /* save old state */
  _FPU_GETCW(oldcw);
  /* Set FPU flags to use double, not double extended,
     with rounding to nearest */
  cw = (_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
  _FPU_SETCW(cw);
  return (unsigned long long) oldcw;
  #else
  return 0;
  #endif /* defined(_MSC_VER) || defined(__MINGW32__) */
#else /* */
  return 0;
#endif /* MW_IS_X86 */
}



/* Print a format string followed by an error code with a string description of the error.
   Somewhere between standard perror() and warn(). Includes user stuff but without the noise of process name
*/
void mwPerror(const char* fmt, ...)
{
    va_list argPtr;

    va_start(argPtr, fmt);
    vfprintf(stderr, fmt, argPtr);
    va_end(argPtr);

    fprintf(stderr, " (%d): %s\n", errno, strerror(errno));
}

#ifdef _WIN32

#if 1
#define ERROR_MSG_LANG MAKELANGID(LANG_ENGLISH, SUBLANG_ENGLISH_US)
#else
/* Default language */
#define ERROR_MSG_LANG MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT)
#endif

/* Like mwPerror, but for Win32 API functions which use GetLastError() */
void mwPerrorW32(const char* fmt, ...)
{
    va_list argPtr;
    LPVOID msgBuf;
    DWORD rc;

    static const DWORD flags = FORMAT_MESSAGE_ALLOCATE_BUFFER
                             | FORMAT_MESSAGE_FROM_SYSTEM
                             | FORMAT_MESSAGE_IGNORE_INSERTS;

    rc = FormatMessage(flags,
                       NULL,
                       GetLastError(),
                       ERROR_MSG_LANG,
                       (LPTSTR) &msgBuf,
                       0,
                       NULL);

    va_start(argPtr, fmt);
    vfprintf(stderr, fmt, argPtr);
    va_end(argPtr);

    fprintf(stderr, " (%ld): %s", GetLastError(), msgBuf);
    LocalFree(msgBuf);
}
#endif /* _WIN32 */


