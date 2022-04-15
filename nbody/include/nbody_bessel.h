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

real_0 besselJ0(const real_0 x);
real_0 besselJ1(const real_0 x);
real_0 besselJ2(const real_0 x);
real_0 besselI0(const real_0 x);
real_0 besselI1(const real_0 x);
real_0 besselK0(const real_0 x);
real_0 besselK1(const real_0 x);
real_0 besselK2(const real_0 x);

real mw_besselJ0(real* a);
real mw_besselJ1(real* a);
real mw_besselJ2(real* a);
real mw_besselI0(real* a);
real mw_besselI1(real* a);
real mw_besselK0(real* a);
real mw_besselK1(real* a);
real mw_besselK2(real* a);

real_0 besselJ0_zero(const int n);
real_0 besselJ1_zero(const int n);
real aExp(const real* k, const real* R, const real* Rd);
real bExp(const real* k, const real* R, const real* Rd);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_BESSEL_H_ */

