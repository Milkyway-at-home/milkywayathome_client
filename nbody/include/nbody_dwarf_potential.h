/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016-2018 Siddhartha Shelton

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

#ifndef _NBODY_DWARF_POTENTIAL_H_
#define _NBODY_DWARF_POTENTIAL_H_

#include <lua.h>
#include "nbody_types.h"
#include "milkyway_math.h"
#include "nbody_potential_types.h"

#ifdef __cplusplus
extern "C" {
#endif

real_0 get_potential(const Dwarf* args, real_0 r);
real_0 get_density(const Dwarf* args, real_0 r);

real_0 get_first_derv_potential(const Dwarf* args, real_0 r);
real_0 get_first_derv_density(const Dwarf* args, real_0 r);
real_0 get_second_derv_potential(const Dwarf* args, real_0 r);
real_0 get_second_derv_density(const Dwarf* args, real_0 r);

real get_potential_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight);
real get_density_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight);

real get_first_derv_potential_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight);
real get_first_derv_density_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight);
real get_second_derv_potential_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight);
real get_second_derv_density_real(const Dwarf* model_light, const Dwarf* model_dark, real* r, mwbool isLight);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_DWARF_POTENTIAL_H_ */

