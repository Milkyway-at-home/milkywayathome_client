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

#ifndef _NBODY_SHOW_H_
#define _NBODY_SHOW_H_

#include "nbody_types.h"
#include "milkyway_show.h"

/* Types -> String */
const char* showBool(mwbool);
const char* showCriterionT(criterion_t);
const char* showSphericalT(spherical_t);
const char* showDiskT(disk_t);
const char* showHaloT(halo_t);
const char* showNBodyStatus(NBodyStatus);

char* showSpherical(const Spherical*);
char* showHalo(const Halo*);
void printHalo(const Halo*);
char* showDisk(const Disk*);
void printDisk(const Disk*);
char* showPotential(const Potential*);
void printPotential(const Potential*);

char* showNBodyCtx(const NBodyCtx*);
void printNBodyCtx(const NBodyCtx*);

char* showHistogramParams(const HistogramParams*);
void printHistogramParams(const HistogramParams*);

char* showBody(const Body* p);
void printBody(const Body* p);

void printBodies(const Body* bs, unsigned int n);

char* showNBodyState(const NBodyState*);
void printNBodyState(const NBodyState*);

char* showNBodyTree(const NBodyTree*);
void printNBodyTree(const NBodyTree*);


#endif /* _NBODY_SHOW_H_ */

