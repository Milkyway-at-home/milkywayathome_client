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

#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "io.h"
#include "lua_type_marshal.h"
#include "lua_nbodyctx.h"

#include "milkyway_util.h"

#include "nbody_scriptable.h"

static void passListTest(lua_State* luaSt)
{
/* C array to Lua */

    warn("Pass list test\n");
    real arst[15];

    int table;
    unsigned int i;

    for (i = 0; i < 15; ++i)
        arst[i] = (real) (i + drand48());

    lua_getglobal(luaSt, "takeList");

    lua_createtable(luaSt, 15, 0);
    table = lua_gettop(luaSt);

    NBodyCtx* lotsOCtx;

    lotsOCtx = mwCalloc(15, sizeof(NBodyCtx));

    for (i = 0; i < 15; ++i)
    {
        //lua_pushvalue(luaSt, arst[i]);
        //lua_pushnumber(luaSt, arst[i]);
        lotsOCtx[i].timestep = 1337.0;
        NBodyCtx* tmp = pushNBodyCtx(luaSt);
        int lol = -1;
        warn("last type = %s\n", lua_typename(luaSt, lol));
        *tmp = lotsOCtx[i];
        lua_pushvalue(luaSt, lol);
        lua_rawseti(luaSt, table, i + 1);

        lua_pop(luaSt, 1);
    }

    if (lua_pcall(luaSt, 1, 0, 0))
        warn("Calling takeList failed\n");
}

static void callTestLists(lua_State* luaSt)
{
    /* the function name */
    lua_getglobal(luaSt, "testLists");

    /* the first argument */
    lua_pushnumber(luaSt, 9.0);

    if (lua_pcall(luaSt, 1, 1, 0))
    {
        warn("Calling testLists failed\n");
        return;
    }

    luaL_checktype(luaSt, -1, LUA_TTABLE);

    int i, n;

    n = luaL_getn(luaSt, -1);  /* get size of table */
    warn("Table is %u\n", n);

    int table = lua_gettop(luaSt);
    /* 1st argument must be a table (t) */
    luaL_checktype(luaSt, table, LUA_TTABLE);
    n = luaL_getn(luaSt, table);  /* get size of table */
    warn("Table size is %u\n", n);

    for (i = 1; i <= n; ++i)
    {
        warn("Reading element %u of %u\n", i, n);
        lua_rawgeti(luaSt, table, i);  /* push t[i] */
        luaL_checktype(luaSt, -1, LUA_TNUMBER);
        real readi = lua_tonumber(luaSt, -1);
        warn("Read in[%u] = %f\n", i, readi);
        lua_pop(luaSt, 1);
    }

}


static void callAdd(lua_State* luaSt)
{
    /* the function name */
    lua_getglobal(luaSt, "add");

    /* the first argument */
    lua_pushnumber(luaSt, 41);

    /* the second argument */
    lua_pushnumber(luaSt, 22);

    /* call the function with 2 arguments, return 1 result */
    if (lua_pcall(luaSt, 2, 1, 0))
    {
        warn("Call to add failed\n");
        return;
    }

    /* get the result */
    int sum = (int) lua_tonumber(luaSt, -1);
    lua_pop(luaSt, 1);

    /* print the result */
    warn("The result is %d\n", sum);
}

static void callTestContext(lua_State* luaSt)
{
    NBodyCtx* ctx;

    /* the function name */
    lua_getglobal(luaSt, "testContext");

    if (lua_pcall(luaSt, 0, 1, 0))
    {
        warn("Calling testContext failed\n");
        return;
    }

    /* get the result */
    ctx = checkNBodyCtx(luaSt, -1);
    lua_pop(luaSt, 1);

    /* print the result */
    warn("CONTEXT TEST %f\n", ctx->timestep);
}

int scriptableArst()
{
    lua_State* luaSt;

    luaSt = lua_open();
    if (!luaSt)
    {
        warn("Opening Lua state failed\n");
        return 1;
    }

    printf("initop = %d\n", lua_gettop(luaSt));

    //lu_load(luaSt, reader

    /* Load various Lua libraries */
    //lua_baselibopen(luaSt);


    luaopen_base(luaSt);
    luaopen_table(luaSt);
#if 0
    luaopen_io(luaSt);
    luaopen_string(luaSt);

    luaopen_math(luaSt);
#endif

    printf("arst top = %d\n", lua_gettop(luaSt));

    registerNBodyCtx(luaSt);

    printf("nbody top = %d\n", lua_gettop(luaSt));

    /* load the script */

    // luaL_dostring
    if (luaL_dofile(luaSt, "add.lua") != 0)
        warn("dofile failed\n");

    printf("dofile top = %d\n", lua_gettop(luaSt));
    callAdd(luaSt);
    callTestContext(luaSt);


    callTestLists(luaSt);
    passListTest(luaSt);

    lua_close(luaSt);

    return 0;
}

