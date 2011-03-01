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
#include "nbody_check_params.h"

static int getNBodyCtxFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeContext");
}

static int getModelsFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeModels");
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
        mw_lua_pcall_warn(luaSt, "Error evaluating script");
        lua_close(luaSt);
        luaSt = NULL;
    }

    free(script);
    return luaSt;
}

static int evaluateContext(lua_State* luaSt, NBodyCtx* ctx)
{
    getNBodyCtxFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating NBodyCtx");
        return 1;
    }

    *ctx = *checkNBodyCtx(luaSt, lua_gettop(luaSt));
    lua_pop(luaSt, 1);

    return contextSanityCheck(ctx);
}

static body* evaluateBodies(lua_State* luaSt, unsigned int* n)
{
    getModelsFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating bodies");
        return NULL;
    }

    return readReturnedModels(luaSt, lua_gettop(luaSt), n);
}

static int setupInitialNBodyState(lua_State* luaSt, NBodyCtx* ctx, NBodyState* st)
{
    body* bodies;
    unsigned int n;

    if (evaluateContext(luaSt, ctx))
        return 1;

    bodies = evaluateBodies(luaSt, &n);
    if (!bodies)
        return 1;

    ctx->nbody = n;

    st->tree.rsize = ctx->treeRSize;
    st->tnow = 0.0;
    st->bodytab = bodies;

    return 0;
}

int setupNBody(const char* filename, NBodyCtx* ctx, NBodyState* st)
{
    int rc;
    lua_State* luaSt;

    luaSt = openNBodyLuaState(filename);
    if (!luaSt)
        return 1;

    rc = setupInitialNBodyState(luaSt, ctx, st);
    lua_close(luaSt);

    return rc;
}


