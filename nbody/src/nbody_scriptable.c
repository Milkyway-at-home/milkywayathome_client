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
#include "lua_body.h"
#include "lua_body_array.h"
#include "nbody_scriptable.h"

#include "milkyway_util.h"


#define TOP_TYPE(st, msg) warn("%s: %s\n", msg, luaL_typename(st, -1));


static void pushRealArray(lua_State* luaSt, const real* arr, int n)
{
    int i, table;

    lua_createtable(luaSt, n, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < n; ++i)
    {
        lua_pushnumber(luaSt, (lua_Number) arr[i]);
        lua_rawseti(luaSt, table, i + 1);
        lua_pop(luaSt, 1);
    }
}

static int pushNBodyCtxArray(lua_State* luaSt, const NBodyCtx* ctxs, int n)
{
    int i, table;

    lua_createtable(luaSt, n, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < n; ++i)
    {
        if (pushNBodyCtx(luaSt, &ctxs[i]))
        {
            warn("Pushing NBodyCtx %d of %d failed\n", i, n);
            lua_pop(luaSt, 1); /* Cleanup table */
            return 1;
        }

        lua_pushvalue(luaSt, -1);
        lua_rawseti(luaSt, table, i + 1);
        lua_pop(luaSt, 1);
    }

    return 1;
}

static real* popRealArray(lua_State* luaSt, int* outN)
{
    real* arr;
    int i, n, table;

    table = lua_gettop(luaSt);
    luaL_checktype(luaSt, table, LUA_TTABLE);
    n = luaL_getn(luaSt, table);  /* get size of table */

    arr = mwMalloc(sizeof(real) * n);
    for (i = 0; i < n; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);  /* push t[i] */
        luaL_checktype(luaSt, -1, LUA_TNUMBER);
        arr[i] = lua_tonumber(luaSt, -1);
        lua_pop(luaSt, 1);
    }

    lua_pop(luaSt, 1);

    if (outN)
        *outN = n;

    return arr;
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

static void callTakeContext(lua_State* luaSt)
{
    NBodyCtx ctx = EMPTY_CTX;

    ctx.timestep = 934.0;
    ctx.orbit_timestep = 12345.0;

    warn("ini = %f, %f\n", ctx.timestep, ctx.orbit_timestep);

    /* the function name */
    lua_getglobal(luaSt, "takeContext");

    if (pushNBodyCtx(luaSt, &ctx))
    {
        warn("Pushing NBodyCtx failed\n");
        return;
    }

    if (lua_pcall(luaSt, 1, 1, 0))
    {
        warn("Calling takeContext failed\n");
        return;
    }

    /* get the result */
    NBodyCtx* ctxResult = checkNBodyCtx(luaSt, -1);
    lua_pop(luaSt, 1);

    /* print the result */
    warn("context modified %f %f\n", ctxResult->timestep, ctxResult->orbit_timestep);
}


static int callProcessContext(lua_State* luaSt, NBodyCtx* ctx)
{
    NBodyCtx* ctxRead;

    lua_getglobal(luaSt, "processContext");
    pushNBodyCtx(luaSt, ctx);

    ctxRead = checkNBodyCtx(luaSt, -1);
    if (!ctxRead)
        return 1;

    *ctx = *ctxRead; /* CHECKME: Who owns this? */

    return 0;
}

static void callTestBodies(lua_State* luaSt)
{
    NBodyLuaBodyArray* bodies = NULL;

    TOP_TYPE(luaSt, "test bodies");

    lua_getglobal(luaSt, "testBodyArray");
    TOP_TYPE(luaSt, "got global");
    lua_call(luaSt, 0, 1);
    TOP_TYPE(luaSt, "call made");
    warn("EVERYTHING IS OK\n");



    /* get the result */
    bodies = checkNBodyLuaBodyArray(luaSt, -1);
    // lua_pop(luaSt, 1);

    /* print the result */
    if (bodies)
        warn("WHEEEEE %u\n", bodies->nBody);
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

    luaopen_base(luaSt);
    luaopen_table(luaSt);
    luaopen_string(luaSt);
    //luaopen_debug(luaSt);
    //luaopen_io(luaSt);

    registerNBodyCtx(luaSt);
    registerNBodyLuaBodyArray(luaSt);

    printf("nbody top = %d\n", lua_gettop(luaSt));



    // luaL_dostring
    if (luaL_dofile(luaSt, "add.lua") != 0)
        warn("dofile failed\n");

    callTestBodies(luaSt);

    printf("dofile top = %d\n", lua_gettop(luaSt));
    callAdd(luaSt);
    callTestContext(luaSt);



    //callTestLists(luaSt);

    printf("soup top = %d\n", lua_gettop(luaSt));
    lua_getglobal(luaSt, "testLists");
    lua_pushnumber(luaSt, 9.0);
    lua_call(luaSt, 1, 1);

    printf("pre pop array top = %d\n", lua_gettop(luaSt));
    real* arst = popRealArray(luaSt, NULL);
    free(arst);


    callTakeContext(luaSt);

    printf("final top = %d\n", lua_gettop(luaSt));
    lua_close(luaSt);

    return 0;
}

