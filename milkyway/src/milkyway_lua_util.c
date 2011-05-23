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
#include <lualib.h>

#include "milkyway_lua_marshal.h"
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

static const luaL_Reg normalLibs[] =
{
    { "",              luaopen_base   },
    { LUA_TABLIBNAME,  luaopen_table  },
    { LUA_STRLIBNAME,  luaopen_string },
 // { LUA_MATHLIBNAME, luaopen_math   },  We replace math with bindings to whatever math we're using
    { NULL, NULL}
};

static const luaL_Reg debugOnlyLibs[] =
{
    { LUA_LOADLIBNAME, luaopen_package },
    { LUA_IOLIBNAME,   luaopen_io      },
    { LUA_OSLIBNAME,   luaopen_os      },
    { LUA_DBLIBNAME,   luaopen_debug   },
    { LUA_MATHLIBNAME, luaopen_math    },
    { NULL, NULL}
};

static void mw_luaL_register(lua_State* luaSt, const luaL_Reg* libs)
{
    while (libs->name)
    {
        lua_pushcfunction(luaSt, libs->func);
        lua_pushstring(luaSt, libs->name);
        lua_call(luaSt, 1, 0);
        ++libs;
    }
}

/* Finer control over which standard libraries are opened for sandboxing */
void mw_lua_openlibs(lua_State* luaSt, mwbool debug)
{
    mw_luaL_register(luaSt, normalLibs);
    if (debug)
        mw_luaL_register(luaSt, debugOnlyLibs);
}

static int doWithArgs(lua_State* luaSt, const char** args, unsigned int nArgs)
{
    unsigned int i;

    if (args)
    {
        for (i = 0; i < nArgs; ++i)
            lua_pushstring(luaSt, args[i]);
    }

    if (lua_pcall(luaSt, nArgs, 0, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating script");
        return 1;
    }

    return 0;
}


int dostringWithArgs(lua_State* luaSt,
                     const char* str,
                     const char** args,
                     unsigned int nArgs)
{
    if (luaL_loadstring(luaSt, str))
        return 1;

    return doWithArgs(luaSt, args, nArgs);
}

int dofileWithArgs(lua_State* luaSt,
                   const char* filename,
                   const char** args,
                   unsigned int nArgs)
{
    if (luaL_loadfile(luaSt, filename))
        return 1;

    return doWithArgs(luaSt, args, nArgs);
}

int mwBindBOINCStatus(lua_State* luaSt)
{
    int isStandalone = TRUE;

    lua_pushboolean(luaSt, BOINC_APPLICATION);
    lua_setglobal(luaSt, "isBOINCApplication");

    #if BOINC_APPLICATION
    isStandalone = boinc_is_standalone();
    #endif

    lua_pushboolean(luaSt, isStandalone);
    lua_setglobal(luaSt, "isStandalone");

    return 0;
}

