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
const char* showBool(const bool);
const char* showCriterionT(const criterion_t);
const char* showSphericalT(const spherical_t);
const char* showDiskT(const disk_t);
const char* showHaloT(const halo_t);
const char* showDwarfModelT(const dwarf_model_t);
char* showSpherical(const Spherical*);
char* showHalo(const Halo*);
char* showDisk(const Disk*);
char* showPotential(const Potential*);
char* showDwarfModel(const DwarfModel*);
char* showInitialConditions(const InitialConditions*);
char* showContext(const NBodyCtx*);
char* showVector(const vector v);
void printContext(const NBodyCtx*);
void printInitialConditions(const InitialConditions*);
void printVector(const vector v);
char* showFitParams(const FitParams* fp);
void printFitParams(const FitParams* fp);

#endif /* _SHOW_H_ */

