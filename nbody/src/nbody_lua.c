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
#include <lauxlib.h>

#include "nbody_types.h"
#include "milkyway_util.h"
#include "orbitintegrator.h"

#include "nbody_lua.h"
#include "nbody_lua_types.h"
#include "nbody_lua_models.h"
#include "lua_type_marshal.h"
#include "nbody_lua_marshal.h"

static int getNBodyCtxFunc(lua_State* luaSt)
{
    return mw_lua_checkglobal(luaSt, "makeContext");
}

static int getPotentialFunc(lua_State* luaSt)
{
    return mw_lua_checkglobal(luaSt, "makePotential");
}

static int getInitialConditionsFunc(lua_State* luaSt)
{
    return mw_lua_checkglobal(luaSt, "makeInitialConditions");
}

static void registerUsedStandardStuff(lua_State* luaSt)
{
    luaopen_base(luaSt);
    luaopen_table(luaSt);
    luaopen_string(luaSt);
    lua_pop(luaSt, 3);
}

static lua_State* openNBodyLuaState(const char* filename)
{
    char* script;
    lua_State* luaSt = NULL;

    script = mwReadFileResolved(filename);
    if (!script)
    {
        perror("Opening Lua script");
        return NULL;
    }

    luaSt = lua_open();
    if (!luaSt)
    {
        warn("Failed to get Lua state\n");
        free(script);
        return NULL;
    }

    registerUsedStandardStuff(luaSt);
    registerNBodyTypes(luaSt);
    registerOtherTypes(luaSt);
    registerPredefinedModelGenerators(luaSt);
    registerModelUtilityFunctions(luaSt);

    if (luaL_dostring(luaSt, script))
    {
        /* TODO: Get error */
        warn("dostring failed\n");
        lua_close(luaSt);
        luaSt = NULL;
    }

    free(script);
    return luaSt;
}

mwbool setupNBody(const char* filename, NBodyCtx* ctx, NBodyState* st, HistogramParams* histParams)
{
    lua_State* luaSt;

    luaSt = openNBodyLuaState(filename);
    if (!luaSt)
        return TRUE;

    getNBodyCtxFunc(luaSt);
    getPotentialFunc(luaSt);
    getInitialConditionsFunc(luaSt);

    mw_lua_checkglobal(luaSt, "arst");
    lua_call(luaSt, 0, 1);
    readReturnedModels(luaSt, lua_gettop(luaSt));

    lua_close(luaSt);

    return FALSE;
}


