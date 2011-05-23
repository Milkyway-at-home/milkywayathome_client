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

#if !defined(_MILKYWAY_LUA_H_INSIDE_) && !defined(MILKYWAY_LUA_COMPILATION)
  #error "Only milkyway_lua.h can be included directly."
#endif

#ifndef _MILKYWAY_LUA_UTIL_H_
#define _MILKYWAY_LUA_UTIL_H_

#include <lua.h>
#include "milkyway_util.h"

int registerUtilityFunctions(lua_State* luaSt);
void mw_lua_openlibs(lua_State *L, mwbool debug);

int dostringWithArgs(lua_State* luaSt, const char* str, const char** args, unsigned int nArgs);
int dofileWithArgs(lua_State* luaSt, const char* filename, const char** args, unsigned int nArgs);

int mwBindBOINCStatus(lua_State* luaSt);

#endif /* _MILKYWAY_LUA_UTIL_H_ */

