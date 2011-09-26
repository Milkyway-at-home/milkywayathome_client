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

#ifndef _MW_CPUID_H_
#define _MW_CPUID_H_
/* returns eax, ebx, ecx, edx in abcd[], using a = eax, c = ecx as input */
void mw_cpuid(int abcd[4], int a, int c);

/* Takes populated array from mw_cpuid() and tests for different instructions */
int mwHasSSE41(const int abcd[4]);
int mwHasSSE3(const int abcd[4]);
int mwHasSSE2(const int abcd[4]);
int mwHasAVX(const int abcd[4]);

#endif /* _MW_CPUID_H_ */


