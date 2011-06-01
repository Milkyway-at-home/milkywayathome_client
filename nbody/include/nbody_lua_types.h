/*
Copyright (C) 2011  Matthew Arsenault

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

#ifndef _NBODY_LUA_TYPES_H_
#define _NBODY_LUA_TYPES_H_

#define _NBODY_LUA_TYPES_H_INSIDE_

#ifdef __cplusplus
extern "C" {
#endif

#include "milkyway_lua.h"

#include "nbody_lua_type_marshal.h"

#include "nbody_lua_nbodyctx.h"
#include "nbody_lua_nbodystate.h"
#include "nbody_lua_body.h"
#include "nbody_lua_halo.h"
#include "nbody_lua_disk.h"
#include "nbody_lua_spherical.h"
#include "nbody_lua_potential.h"
#include "nbody_lua_histogram_params.h"


#include <lua.h>

void registerNBodyTypes(lua_State* luaSt);
NBodyStatus readNBodyStatus(lua_State* luaSt, const char* name);

#ifdef __cplusplus
}
#endif

#undef _NBODY_LUA_TYPES_H_INSIDE_

#endif /* _NBODY_LUA_TYPES_H_ */

