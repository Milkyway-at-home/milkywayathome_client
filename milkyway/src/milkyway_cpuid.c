/*
 *  Copyright (c) 2011 Matthew Arsenault
 *  Copyright (c) 2011 Rensselaer Polytechnic Institute
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

#include "milkyway_cpuid.h"
#include "milkyway_util.h"

#if defined(__APPLE__) && MW_IS_X86
  #include <sys/param.h>
  #include <sys/sysctl.h>
#endif

#define bit_CMPXCHG8B (1 << 8)
#define bit_CMOV (1 << 15)
#define bit_MMX (1 << 23)
#define bit_SSE (1 << 25)
#define bit_SSE2 (1 << 26)
#define bit_SSE3 (1 << 0)
#define bit_SSE41 (1 << 19)
#define bit_AVX (1 << 28)
#define bit_CMPXCHG16B (1 << 13)
#define bit_3DNOW (1 << 31)
#define bit_3DNOWP (1 << 30)
#define bit_LM (1 << 29)


#if MW_IS_X86

#ifndef _WIN32
void mw_cpuid(int abcd[4], int a, int c)
{
    abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;

  #if defined(__i386__) && defined(__PIC__)
    __asm__ volatile(
        "pushl %%ebx     \n\t" /* save %ebx to prevent break with PIC */
        "cpuid           \n\t" /* input eax, ecx. output eax, ebx, ecx, edx */
        "movl %%ebx, %1  \n\t" /* save what cpuid just put in %ebx */
        "popl %%ebx      \n\t" /* restore the old %ebx */
        : "=a" (abcd[0]), "=r" (abcd[1]), "=c" (abcd[2]), "=d" (abcd[3])
        : "a" (a), "c" (c)
        : "cc"
        );
  #else
    __asm__ volatile(
        "cpuid           \n\t" /* input eax, ecx. output eax, ebx, ecx, edx */
        : "=a" (abcd[0]), "=b" (abcd[1]), "=c" (abcd[2]), "=d" (abcd[3])
        : "a" (a), "c" (c)
        );
  #endif
}

#else /* _WIN32 */

void mw_cpuid(int abcd[4], int a, int c)
{
    abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
    __cpuid(abcd, 0);
    if (abcd[0] >= 1) /* Is this really necessary? */
    {
        __cpuid(abcd, a);
    }
    else
    {
        abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
    }
}
#endif /* _WIN32 */

#else /* !MW_IS_X86 */

void mw_cpuid(int abcd[4], int a, int c)
{
    (void) a, (void) c;

    abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
}

#endif /* MW_IS_X86 */

int mwHasAVX(const int abcd[4])
{
    return !!(abcd[2] & bit_AVX);
}

int mwHasSSE41(const int abcd[4])
{
    return !!(abcd[2] & bit_SSE41);
}

int mwHasSSE3(const int abcd[4])
{
    return !!(abcd[2] & bit_SSE3);
}

int mwHasSSE2(const int abcd[4])
{
    return !!(abcd[3] & bit_SSE2);
}



#if defined(_WIN32)
int mwOSHasAVXSupport(void)
{
    OSVERSIONINFOEX vInfo;

    vInfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
    if (!GetVersionEx(&vInfo))
    {
        mwPerrorW32("Error getting Windows version info");
        return FALSE;
    }

    /* Windows 7 SP1 or greater required. Can't find a real way to check. */
    return (vInfo.dwMajorVersion > 6)
        || (vInfo.dwMajorVersion == 6 && vInfo.dwMinorVersion >= 1 && vInfo.wServicePackMajor >= 1);
}

#elif defined(__linux__)

/* Scan through /proc/cpuinfo to see if the kernel supports AVX */
int mwOSHasAVXSupport(void)
{
    FILE* f;
    char lineBuf[4096];
    int support = FALSE;

    f = fopen("/proc/cpuinfo", "r");
    if (!f)
    {
        return FALSE;
    }

    while (fgets(lineBuf, (int) sizeof(lineBuf), f))
    {
        if (strncmp(lineBuf, "flags", sizeof(lineBuf)))
        {
            if (strstr(lineBuf, " avx "))
            {
                support = TRUE;
                break;
            }
        }
    }

    fclose(f);

    return support;
}

#elif defined(__APPLE__) && MW_IS_X86

int mwOSHasAVXSupport(void)
{
    size_t len;
    char* kernelVersion;
    int major = 0;
    int minor = 0;
    int patchLevel = 0;
    static int mib[2] = { CTL_KERN, KERN_OSRELEASE };

    if (sysctl(mib, 2, NULL, &len, NULL, 0) < 0)
    {
        mwPerror("sysctl kernel version size");
        return FALSE;
    }

    kernelVersion = malloc(len * sizeof(char));
    if (!kernelVersion)
    {
        mwPerror("malloc");
        return FALSE;
    }

    if (sysctl(mib, 2, kernelVersion, &len, NULL, 0) < 0)
    {
        mwPerror("sysctl kernel version");
        free(kernelVersion);
        return FALSE;
    }

    if (sscanf(kernelVersion, "%d.%d.%d", &major, &minor, &patchLevel) != 3)
    {
        major = minor = patchLevel = 0;
        /* Old examples seem to say it was 2 digits */
        if (sscanf(kernelVersion, "%d.%d", &major, &minor) != 2)
        {
            major = minor = patchLevel = 0;
        }
    }

    free(kernelVersion);

   /* AVX support added in 10.6.8, which has kernel 10.8 */
    return (major >= 10 && minor >= 8);
}

#else

int mwOSHasAVXSupport()
{
    return FALSE;
}

#endif /* _WIN32 */



