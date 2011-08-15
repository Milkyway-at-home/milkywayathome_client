/*
Copyright (C) 2011  Matthew Arsenault
Copyright (C) 2011  Rensselaer Polytechnic Institute

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

#include "mw_cpuid.h"
#include "milkyway_util.h"

#define bit_CMPXCHG8B (1 << 8)
#define bit_CMOV (1 << 15)
#define bit_MMX (1 << 23)
#define bit_SSE (1 << 25)
#define bit_SSE2 (1 << 26)
#define bit_SSE3 (1 << 0)
#define bit_SSE41 (1 << 19)
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
    abcd[0] = abcd[1] = abcd[2] = abcd[3] = 0;
}

#endif /* MW_IS_X86 */

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

