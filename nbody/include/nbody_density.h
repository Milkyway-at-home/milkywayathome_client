/*
 *  Copyright (c) 2018-2021 Rensselaer Polytechnic Institute
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
#ifndef _NBODY_DENSITY_H_
#define _NBODY_DENSITY_H_

#include "nbody_types.h"

#ifdef __cplusplus
extern "C" {
#endif

real nbExtDensity(const Potential* pot, mwvector pos, real time);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_DENSITY_H_ */

