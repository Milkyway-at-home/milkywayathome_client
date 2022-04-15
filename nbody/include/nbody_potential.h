/*
 *  Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 *  Copyright (c) 2010-2011 Matthew Arsenault
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

#ifndef _NBODY_POTENTIAL_H_
#define _NBODY_POTENTIAL_H_

#include "nbody_types.h"

#ifdef __cplusplus
extern "C" {
#endif

mwvector nbExtAcceleration(const Potential* pot, mwvector pos, real time);
mwvector pointAccel(const mwvector pos, const mwvector pos1, const real mass);
mwvector plummerAccel(const mwvector pos, const mwvector pos1, const real mass, const real scale);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_POTENTIAL_H_ */

