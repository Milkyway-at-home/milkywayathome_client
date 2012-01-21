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

#ifndef _MILKYWAY_LUA_UTIL_H_
#define _MILKYWAY_LUA_UTIL_H_

#include <lua.h>
#include "milkyway_util.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Print the error from pcall assumed to be the top of the stack, and remove the error */
#define mw_lua_perror(luaSt, msg, ...)                                  \
    do                                                                  \
    {                                                                   \
        fprintf(stderr, msg ": %s \n",                                  \
                ##__VA_ARGS__, lua_tostring(luaSt, -1));                \
        lua_pop(luaSt, 1);                                              \
    }                                                                   \
    while (0)                                                           \


lua_State* mw_lua_newstate(void);

int registerUtilityFunctions(lua_State* luaSt);
void mw_lua_openlibs(lua_State *L, mwbool debug);

int dostringWithArgs(lua_State* luaSt, const char* str, const char** args, unsigned int nArgs);
int dofileWithArgs(lua_State* luaSt, const char* filename, const char** args, unsigned int nArgs);

int mwBindBOINCStatus(lua_State* luaSt);

#ifdef __cplusplus
}
#endif

#endif /* _MILKYWAY_LUA_UTIL_H_ */

