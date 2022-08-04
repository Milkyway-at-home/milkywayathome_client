/*
 * Copyright (c) 2018 Eric Mendelsohn
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
 *
 * The bessel function code was pulled directly from Numerical Recipes - The Art of Scientific
 * Computing: Third Edition, by Willian H. Press, Saul A. Teukolsky,
 * William T. Vetterling, and Brian P. Flannery: pgs.280-281.
 * 
 * The Webnotes for bessel code were pulled from <http://numerical.recipes/webnotes/nr3web7.pdf>.
 */

#ifndef _NBODY_BESSEL_H_
#define _NBODY_BESSEL_H_

#include "nbody_types.h"
#include "nbody.h"
#include "milkyway_util.h"


#ifdef __cplusplus
extern "C" {
#endif

real besselJ0(const real x);
real besselJ1(const real x);
real besselJ2(const real x);
real besselI0(const real x);
real besselI1(const real x);
real besselK0(const real x);
real besselK1(const real x);
real besselK2(const real x);
real besselJ0_zero(const int n);
real besselJ1_zero(const int n);
real aExp(const real k, const real R, const real Rd);
real bExp(const real k, const real R, const real Rd);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_BESSEL_H_ */

