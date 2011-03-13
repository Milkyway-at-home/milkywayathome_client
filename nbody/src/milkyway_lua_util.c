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

#include "milkyway_lua_marshal.h"
#include "nbody_lua_types.h"
#include "milkyway_util.h"
#include "milkyway_lua_util.h"

/* Map over values in a table */
static int luaMap(lua_State* luaSt)
{
    int i, n, f = 1, table = 2, newTable = 3;

    if (lua_gettop(luaSt) != 2)
        luaL_argerror(luaSt, 2, "Expected 2 arguments");

    mw_lua_checktable(luaSt, table);
    mw_lua_checkfunction(luaSt, f);

    n = luaL_getn(luaSt, table);

    lua_newtable(luaSt);
    for (i = 0; i < n; ++i)
    {
        lua_pushvalue(luaSt, f);
        lua_rawgeti(luaSt, table, i + 1);
        lua_call(luaSt, 1, 1);
        lua_rawseti(luaSt, newTable, i + 1);
    }

    return 1;
}

/* FIXME?: can be side effecting on accumulator */
static int luaFoldl(lua_State* luaSt)
{
    int i, n, acc, f = 1, ini = 2, table = 3;

    if (lua_gettop(luaSt) != 3)
        luaL_argerror(luaSt, 3, "Expected 3 arguments");

    mw_lua_checkfunction(luaSt, f);
    /* accumulator can be anything */
    mw_lua_checktable(luaSt, table);

    n = luaL_getn(luaSt, table);

    lua_pushvalue(luaSt, ini);
    acc = lua_gettop(luaSt);

    for (i = 0; i < n; ++i)
    {
        lua_pushvalue(luaSt, f);
        lua_pushvalue(luaSt, acc);
        lua_rawgeti(luaSt, table, i + 1);
        lua_call(luaSt, 2, 1);
        lua_replace(luaSt, acc);
    }

    return 1;
}

static int luaZipWith(lua_State* luaSt)
{
    int i, n, f = 1, tableA = 2, tableB = 3, newTable = 4;

    if (lua_gettop(luaSt) != 3)
        luaL_argerror(luaSt, 3, "Expected 3 arguments");

    mw_lua_checkfunction(luaSt, f);
    mw_lua_checktable(luaSt, tableA);
    mw_lua_checktable(luaSt, tableB);

    n = MIN(luaL_getn(luaSt, tableA), luaL_getn(luaSt, tableB));

    lua_newtable(luaSt);
    for (i = 0; i < n; ++i)
    {
        lua_pushvalue(luaSt, f);
        lua_rawgeti(luaSt, tableA, i + 1);
        lua_rawgeti(luaSt, tableB, i + 1);
        lua_call(luaSt, 2, 1);
        lua_rawseti(luaSt, newTable, i + 1);
    }

    return 1;
}

int registerUtilityFunctions(lua_State* luaSt)
{
    lua_register(luaSt, "map", luaMap);
    lua_register(luaSt, "foldl", luaFoldl);
    lua_register(luaSt, "zipWith", luaZipWith);

    return 0;
}

