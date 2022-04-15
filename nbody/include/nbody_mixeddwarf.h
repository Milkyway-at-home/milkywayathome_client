/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright 2016 Siddhartha Shelton

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

#ifndef _NBODY_MIXEDDWARF_H_
#define _NBODY_MIXEDDWARF_H_

#include <lua.h>
#include "nbody_types.h"
#include "nbody_potential_types.h"

#ifdef __cplusplus
extern "C" {
#endif

int nbGenerateMixedDwarfCore_TESTVER(mwvector* pos, mwvector* vel, real* bodyMasses, dsfmt_t* prng, unsigned int nbody, 
                                     Dwarf* comp1,  Dwarf* comp2, mwvector rShift, mwvector vShift);

int nbGenerateMixedDwarf(lua_State* luaSt);
void registerGenerateMixedDwarf(lua_State* luaSt);

#ifdef __cplusplus
}
#endif

#endif /* _NBODY_MIXEDDWARF_H_ */

