/*
 *  Copyright (c) 2011 Matthew Arsenault
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

#if !defined(_MILKYWAY_LUA_H_INSIDE_) && !defined(MILKYWAY_LUA_COMPILATION)
  #error "Only milkyway_lua.h can be included directly."
#endif

#ifndef _MILKYWAY_LUA_DSFMT_H_
#define _MILKYWAY_LUA_DSFMT_H_

#include <lua.h>
#include <dSFMT.h>

#define DSFMT_TYPE "DSFMT"


dsfmt_t* checkDSFMT(lua_State* luaSt, int idx);
int pushDSFMT(lua_State* luaSt, const dsfmt_t* state);
int registerDSFMT(lua_State* luaSt);

#endif /* _MILKYWAY_LUA_DSFMT_H_ */

