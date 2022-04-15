/*
 * Copyright (c) 2020 Eric Mendelsohn
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

#ifndef _NBODY_VIRIAL_H_
#define _NBODY_VIRIAL_H_

#include "nbody_types.h"
#include "nbody.h"
#include "milkyway_util.h"


#ifdef __cplusplus
extern "C" {
#endif

real nbCalculateVirial(real a_b, real a_d, real M_b, real M_d);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_VIRIAL_H_ */
