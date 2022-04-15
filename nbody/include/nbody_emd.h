/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NBODY_EMD_H_
#define _NBODY_EMD_H_

#include "milkyway_extra.h"
#include "nbody_types.h"


typedef struct
{
    real_0 weight; /* Normalized Weight */
    real_0 lambda;    /* Lambda Position */
    real_0 beta; /* Beta Position */
} WeightPos;



#ifdef __cplusplus
extern "C" {
#endif

real_0 emdCalc(const real_0* RESTRICT signature_arr1,
              const real_0* RESTRICT signature_arr2,
              unsigned int size1,
              unsigned int size2,
              real_0* RESTRICT lower_bound);

real nbMatchEMD(const MainStruct* data, const MainStruct* histogram);

real_0 nbWorstCaseEMD(const NBodyHistogram* hist  );

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_EMD_H_ */

