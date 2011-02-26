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
#include "nbody_scriptable.h"
#include "nbody_lua_types.h"
#include "nbody_lua_functions.h"

#include "milkyway_util.h"
#include "show.h"


#define TOP_TYPE(st, msg) warn("%s: %s\n", msg, luaL_typename(st, -1));
#define WHEREAMI(msg, st) warn("%s: %d\n", msg, lua_gettop(luaSt));

#if 0
int l_map(lua_State* L)
{
    int i, n;

    /* 1st argument must be a table (t) */
    luaL_checktype(L, 1, LUA_TTABLE);

    /* 2nd argument must be a function (f) */
    luaL_checktype(L, 2, LUA_TFUNCTION);

    n = luaL_getn(L, 1);  /* get size of table */
    for (i = 1; i <= n; ++i)
    {
        lua_pushvalue(L, 2);   /* push f */
        lua_rawgeti(L, 1, i);  /* push t[i] */
        lua_call(L, 1, 1);     /* call f(t[i]) */
        lua_rawseti(L, 1, i);  /* t[i] = result */
    }

    return 0;  /* no results */
}
#endif


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
        pushNBodyCtx(luaSt, &ctxs[i]);
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

    lua_call(luaSt, 0, 1);

    #if 0
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        warn("Calling testContext failed\n");
        return;
    }
    #endif

    /* get the result */
    ctx = checkNBodyCtx(luaSt, -1);
    printNBodyCtx(ctx);
    lua_pop(luaSt, 1);
}

static void callTestInitialConditions(lua_State* luaSt)
{
    InitialConditions* ic;

    lua_getglobal(luaSt, "testInitialConditions");
    lua_call(luaSt, 0, 1);

    ic = checkInitialConditions(luaSt, -1);
    printInitialConditions(ic);
    lua_pop(luaSt, 1);
}

static void callTestGeneratePlummer(lua_State* luaSt)
{
    //TOP_TYPE(luaSt, "Generate plummer call")
    warn("Generate plummer call: %d\n", lua_gettop(luaSt));
    lua_getglobal(luaSt, "testGeneratePlummer");
    lua_call(luaSt, 0, 0);
}

static void callTestHalo(lua_State* luaSt)
{
    Halo* h = NULL;

    lua_getglobal(luaSt, "testHalo");
    lua_pushliteral(luaSt, "nfw");
    lua_call(luaSt, 1, 1);

    h = checkHalo(luaSt, -1);
    printHalo(h);
    lua_pop(luaSt, 1);
}

static void callTestDisk(lua_State* luaSt)
{
    Disk* d = NULL;

    lua_getglobal(luaSt, "testDisk");
    lua_pushliteral(luaSt, "exponential");
    lua_call(luaSt, 1, 1);

    d = checkDisk(luaSt, -1);
    printDisk(d);
    lua_pop(luaSt, 1);
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
    NBodyCtx* ctxResult;

    ctxResult = checkNBodyCtx(luaSt, -1);
    /* print the result */
    warn("context modified %f %f\n", ctxResult->timestep, ctxResult->orbit_timestep);

    lua_pop(luaSt, 1);

}


static int callProcessContext(lua_State* luaSt, NBodyCtx* ctx)
{
    NBodyCtx* ctxRead;

    lua_getglobal(luaSt, "processContext");
    pushNBodyCtx(luaSt, ctx);

    ctxRead = checkNBodyCtx(luaSt, -1);
    *ctx = *ctxRead; /* CHECKME: Who owns this? */
    lua_pop(luaSt, 1);

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
    //bodies = checkNBodyLuaBodyArray(luaSt, -1);
    // lua_pop(luaSt, 1);

    body* b = checkBody(luaSt, -1);
    printBody(b);
    lua_pop(luaSt, 1);


    /* print the result */
    if (bodies)
        warn("WHEEEEE %u\n", bodies->nBody);
}

static void callTestVector(lua_State* luaSt)
{
    lua_getglobal(luaSt, "testVector");
    lua_call(luaSt, 0, 1);

    mwvector arstv = *checkVector(luaSt, -1);
    warn("arst = %f, %f, %f, %f\n", arstv.x, arstv.y, arstv.z, arstv.w);
    lua_pop(luaSt, 1);
}

static void callTestDSFMT(lua_State* luaSt)
{
    dsfmt_t* rSt = NULL;

    TOP_TYPE(luaSt, "test dsfmt");

    lua_getglobal(luaSt, "testDSFMT");
    lua_call(luaSt, 0, 1);

    rSt = checkDSFMT(luaSt, -1);
    lua_pop(luaSt, 1);
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
    lua_pop(luaSt, 3);
    //luaopen_debug(luaSt);
    //luaopen_io(luaSt);
    WHEREAMI("opened standard libraries", luaSt);

    TOP_TYPE(luaSt, "opened std libs");

    registerVector(luaSt);
    registerBody(luaSt);
    registerHalo(luaSt);
    WHEREAMI("some registered", luaSt);
    registerDisk(luaSt);
    registerInitialConditions(luaSt);
    registerNBodyLuaBodyArray(luaSt);
    registerNBodyCtx(luaSt);
    registerDSFMT(luaSt);
    registerPredefinedModelGenerators(luaSt);

    //lua_pop(luaSt, 11);

    WHEREAMI("Everything registered", luaSt);

    // luaL_dostring
    if (luaL_dofile(luaSt, "add.lua") != 0)
        warn("dofile failed\n");

    callTestHalo(luaSt);
    WHEREAMI("callTestHalo", luaSt);
    callTestDisk(luaSt);
    WHEREAMI("callTestDisk", luaSt);
    callTestBodies(luaSt);
    WHEREAMI("callTestBodies", luaSt);

    callTestGeneratePlummer(luaSt);
    WHEREAMI("callTestGeneratePlummer A", luaSt);


    WHEREAMI("callTestDSFMT", luaSt);
    callTestDSFMT(luaSt);

    WHEREAMI("callAdd", luaSt);
    callAdd(luaSt);

    WHEREAMI("callTestVector", luaSt);
    callTestVector(luaSt);

    WHEREAMI("callTestContext", luaSt);
    callTestContext(luaSt);

    WHEREAMI("callTestInitialConditions", luaSt);
    callTestInitialConditions(luaSt);

    WHEREAMI("callTestGeneratePlummer", luaSt);
    callTestGeneratePlummer(luaSt);




    WHEREAMI("call test lists", luaSt);
    callTestLists(luaSt);

    WHEREAMI("pre pop array top", luaSt);
    real* arst = popRealArray(luaSt, NULL);
    free(arst);

    WHEREAMI("callTakeContext", luaSt);
    callTakeContext(luaSt);


    WHEREAMI("Final", luaSt);
    lua_close(luaSt);

    return 0;
}

