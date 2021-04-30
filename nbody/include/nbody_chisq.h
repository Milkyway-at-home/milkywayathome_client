/*
 * Copyright (c) 2010, 2011 Ben Willett
 * Copyright (c) 2010, 2011 Matthew Arsenault
 * Copyright (c) 2010, 2011 Rensselaer Polytechnic Institute.
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

#ifndef _NBODY_CHISQ_H_
#define _NBODY_CHISQ_H_

#include "nbody_types.h"
#include "nbody.h"


#ifdef __cplusplus
extern "C" {
#endif

real nbCalcChisq(const NBodyHistogram* data,        /* Data histogram */
                   const NBodyHistogram* histogram,   /* Generated histogram */
                   NBodyLikelihoodMethod method);


#ifdef __cplusplus
}
#endif

#endif /* _NBODY_CHISQ_H_ */

