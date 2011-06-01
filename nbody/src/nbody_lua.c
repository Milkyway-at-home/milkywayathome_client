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

#include "nbody_lua.h"
#include "nbody_lua_types.h"
#include "nbody_lua_models.h"
#include "milkyway_lua.h"
#include "nbody_check_params.h"

static int getNBodyCtxFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeContext");
}

static int getNBodyPotentialFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makePotential");
}

static int getHistogramFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeHistogram");
}

static int getBodiesFunc(lua_State* luaSt)
{
    return mw_lua_checkglobalfunction(luaSt, "makeBodies");
}

static int bindArgSeed(lua_State* luaSt, const NBodyFlags* nbf)
{
    if (nbf->setSeed)
    {
        lua_pushinteger(luaSt, nbf->setSeed);
        lua_setglobal(luaSt, "argSeed");
    }

    return 0;
}

/* Open a lua_State and load the stuff we define, but do not run anything */
lua_State* nbodyLuaOpen(mwbool debug)
{
    lua_State* luaSt;

    luaSt = lua_open();
    if (!luaSt)
    {
        warn("Failed to get Lua state\n");
        return NULL;
    }

    mw_lua_openlibs(luaSt, debug);

    registerNBodyTypes(luaSt);
    mwRegisterTypes(luaSt);

    registerPredefinedModelGenerators(luaSt);
    registerSphericalKinds(luaSt);
    registerDiskKinds(luaSt);
    registerHaloKinds(luaSt);

    registerModelUtilityFunctions(luaSt);
    registerUtilityFunctions(luaSt);

    return luaSt;
}

/* Open a lua_State, bind run information such as server arguments and
 * BOINC status, and evaluate input script. */
static lua_State* nbodyOpenLuaStateWithScript(const NBodyFlags* nbf)
{
    char* script;
    lua_State* luaSt;

    luaSt = nbodyLuaOpen(nbf->debugLuaLibs);
    if (!luaSt)
        return NULL;

    bindArgSeed(luaSt, nbf);
    mwBindBOINCStatus(luaSt);

    script = mwReadFileResolved(nbf->inputFile);
    if (!script)
    {
        perror("Opening Lua script");
        lua_close(luaSt);
        return NULL;
    }

    if (dostringWithArgs(luaSt, script, nbf->forwardedArgs, nbf->numForwardedArgs))
    {
        mw_lua_pcall_warn(luaSt, "Error loading Lua script '%s'", nbf->inputFile);
        lua_close(luaSt);
        luaSt = NULL;
    }

    free(script);
    return luaSt;
}

static int evaluateContext(lua_State* luaSt, NBodyCtx* ctx)
{
    NBodyCtx* tmp;

    getNBodyCtxFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating NBodyCtx");
        return 1;
    }

    tmp = expectNBodyCtx(luaSt, lua_gettop(luaSt));
    if (!tmp)
    {
        lua_pop(luaSt, 1);
        return 1;
    }

    *ctx = *tmp;
    lua_pop(luaSt, 1);

    return checkNBodyCtxConstants(ctx);
}

static int evaluatePotential(lua_State* luaSt, NBodyCtx* ctx)
{
    int top;
    Potential* tmp;

    getNBodyPotentialFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating Potential");
        return 1;
    }

    top = lua_gettop(luaSt);
    if (lua_isnoneornil(luaSt, top))
    {
        ctx->potentialType = EXTERNAL_POTENTIAL_NONE;
        lua_pop(luaSt, 1);
        return FALSE;
    }
    else if (lua_isfunction(luaSt, top))
    {
        mw_panic("Feature not implemented\n");
    }
    else
    {
        tmp = expectPotential(luaSt, top);
        if (!tmp)
        {
            lua_pop(luaSt, 1);
            return 1;
        }

        ctx->potentialType = EXTERNAL_POTENTIAL_DEFAULT;
        ctx->pot = *tmp;

        lua_pop(luaSt, 1);
        return checkPotentialConstants(&ctx->pot);
    }

    assert(0);
    return 0;
}

static int evaluateHistogram(lua_State* luaSt, HistogramParams* hp)
{
    HistogramParams* tmp;
    int rc = 0;

    getHistogramFunc(luaSt);
    if (lua_pcall(luaSt, 0, 1, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating HistogramParams");
        return 1;
    }

    tmp = expectHistogramParams(luaSt, lua_gettop(luaSt));
    if (!tmp)
        rc = 1;
    else
        *hp = *tmp;
    lua_pop(luaSt, 1);
    return rc;
}

static Body* evaluateBodies(lua_State* luaSt, const NBodyCtx* ctx, unsigned int* n)
{
    int level, nResults;

    level = lua_gettop(luaSt);

    getBodiesFunc(luaSt);
    pushNBodyCtx(luaSt, ctx);

    if (ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT)
        pushPotential(luaSt, &ctx->pot);
    else
        lua_pushnil(luaSt);

    if (lua_pcall(luaSt, 2, LUA_MULTRET, 0))
    {
        mw_lua_pcall_warn(luaSt, "Error evaluating bodies");
        return NULL;
    }

    nResults = lua_gettop(luaSt) - level;

    return readModels(luaSt, nResults, n);
}

static int evaluateInitialNBodyState(lua_State* luaSt, NBodyCtx* ctx, NBodyState* st, HistogramParams* hp)
{
    Body* bodies;
    unsigned int nbody;

    if (evaluateContext(luaSt, ctx))
        return 1;

    if (evaluatePotential(luaSt, ctx))
        return 1;

    if (evaluateHistogram(luaSt, hp))
        return 1;

    bodies = evaluateBodies(luaSt, ctx, &nbody);
    if (!bodies)
        return 1;

    setInitialNBodyState(st, ctx, bodies, nbody);

    return 0;
}

int setupNBody(NBodyCtx* ctx, NBodyState* st, HistogramParams* hp, const NBodyFlags* nbf)
{
    int rc;
    lua_State* luaSt;

    luaSt = nbodyOpenLuaStateWithScript(nbf);
    if (!luaSt)
        return 1;

    rc = evaluateInitialNBodyState(luaSt, ctx, st, hp);
    lua_close(luaSt);

    return rc;
}


