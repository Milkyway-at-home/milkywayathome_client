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

#ifndef _MILKYWAY_CPUID_H_
#define _MILKYWAY_CPUID_H_

#ifdef __cplusplus
extern "C" {
#endif

/* returns eax, ebx, ecx, edx in abcd[], using a = eax, c = ecx as input */
void mw_cpuid(int abcd[4], int a, int c);

/* Takes populated array from mw_cpuid() and tests for different instructions */
int mwHasSSE41(const int abcd[4]);
int mwHasSSE3(const int abcd[4]);
int mwHasSSE2(const int abcd[4]);
int mwHasAVX(const int abcd[4]);

int mwOSHasAVXSupport(void);

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_CPUID_H_ */


