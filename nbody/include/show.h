/* Copyright 2010 Matthew Arsenault, Travis Desell, Dave Przybylo,
Nathan Cole, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#ifndef _SHOW_H_
#define _SHOW_H_

#define _GNU_SOURCE

#include "nbody_types.h"

/* Types -> String */
const char* showBool(const mwbool);
const char* showCriterionT(const criterion_t);
const char* showSphericalT(const spherical_t);
const char* showDiskT(const disk_t);
const char* showHaloT(const halo_t);
char* showSpherical(const Spherical*);
char* showHalo(const Halo*);
void printHalo(const Halo*);
char* showDisk(const Disk*);
void printDisk(const Disk*);
char* showPotential(const Potential*);
void printPotential(const Potential*);

char* showInitialConditions(const InitialConditions*);
char* showNBodyCtx(const NBodyCtx*);
void printNBodyCtx(const NBodyCtx*);
char* showVector(const mwvector);

void printInitialConditions(const InitialConditions*);
void printVector(const mwvector v);

char* showHistogramParams(const HistogramParams*);
void printHistogramParams(const HistogramParams*);

char* showBody(const bodyptr p);
void printBody(const bodyptr p);

#endif /* _SHOW_H_ */

