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

#include <lua.h>
#include <lauxlib.h>

#include "lua_type_marshal.h"
#include "nbody_lua_types.h"
#include "nbody_types.h"
#include "milkyway_util.h"
#include "nbody_lua_marshal.h"


body* popBodyArray(lua_State* luaSt, int table, int* outN)
{
    body* arr;
    int i, n;

    warn("popping\n");

    luaL_checktype(luaSt, table, LUA_TTABLE);
    n = luaL_getn(luaSt, table);  /* get size of table */

    warn("arsting\n");
    double t1 = mwGetTime();

    arr = (body*) mwMalloc(sizeof(body) * n);
    for (i = 0; i < n; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);  /* push t[i] */
        arr[i] = *checkBody(luaSt, -1);
        //printBody(&arr[i]);
        lua_pop(luaSt, 1);
    }

    double t2 = mwGetTime();

    warn("time to pop = %f\n", t2 - t1);

    if (outN)
        *outN = n;

    return arr;
}

void readReturnedModels()
{

}


