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

#include "nbody_types.h"
#include "milkyway_util.h"
#include "orbitintegrator.h"

#include "nbody_lua_functions.h"
#include "nbody_lua_types.h"
#include "lua_type_marshal.h"
#include "plummer.h"
#include "nbody_lua_marshal.h"


void registerPredefinedModelGenerators(lua_State* luaSt)
{
    registerGeneratePlummer(luaSt);
}

/* For using a combination of light and dark models to generate timestep */
static real plummerTimestepIntegral(real smalla, real biga, real Md, real step)
{
    /* Calculate the enclosed mass of the big sphere within the little sphere's scale length */
    real encMass, val, r;

    encMass = 0.0;
    for (r = 0.0; r <= smalla; r += step)
    {
        val = sqr(r) / mw_pow(sqr(r) + sqr(biga), 2.5);
        encMass += val * step;
    }
    encMass *= 3.0 * Md * sqr(biga);

    return encMass;
}

static int luaPlummerTimestepIntegral(lua_State* luaSt)
{
    int nArgs;
    static real smalla = 0.0, biga = 0.0, Md = 0.0, encMass = 0.0;
    static real step = 1.0e-5;

    static const MWNamedArg argTable[] =
        {
            { "smalla", LUA_TNUMBER, NULL, TRUE,  &smalla },
            { "biga",   LUA_TNUMBER, NULL, TRUE,  &biga   },
            { "Md",     LUA_TNUMBER, NULL, TRUE,  &Md     },
            { "step",   LUA_TNUMBER, NULL, FALSE, &step   },
            END_MW_NAMED_ARG
        };

    nArgs = lua_gettop(luaSt);
    if (nArgs == 1 && lua_istable(luaSt, 1))
    {
        handleNamedArgumentTable(luaSt, argTable, 1);
    }
    else if (nArgs == 3 || nArgs == 4)
    {
        smalla = luaL_checknumber(luaSt, 1);
        biga = luaL_checknumber(luaSt, 2);
        Md = luaL_checknumber(luaSt, 3);
        step = luaL_optnumber(luaSt, 4, step);
    }
    else
    {
        return luaL_argerror(luaSt, 1, "Expected 1, 3 or 4 arguments");
    }


    /* Make sure the bounds / step are OK so that this integral will be sure to complete */
    if (mwCheckNormalPosNum(smalla))
        return luaL_argerror(luaSt, 1, "Invalid small radius");
    if (mwCheckNormalPosNum(biga))
        return luaL_argerror(luaSt, 2, "Invalid big radius");
    if (mwCheckNormalPosNum(step))
        return luaL_argerror(luaSt, 4, "Invalid step argument");

    encMass = plummerTimestepIntegral(smalla, biga, Md, step);
    lua_pushnumber(luaSt, encMass);

    return 1;
}

static void registerPlummerTimestepIntegral(lua_State* luaSt)
{
    lua_pushcfunction(luaSt, luaPlummerTimestepIntegral);
    lua_setglobal(luaSt, "plummerTimestepIntegral");
}

static int luaReverseOrbit(lua_State* luaSt)
{
    real dt, tstop;
    const Potential* pot = NULL;
    const InitialConditions* ic = NULL;
    InitialConditions icNew;

    const MWNamedArg argTable[] =
        {
            { "potential",         LUA_TUSERDATA, POTENTIAL_TYPE,          TRUE, &pot   },
            { "initialConditions", LUA_TUSERDATA, INITIAL_CONDITIONS_TYPE, TRUE, &ic    },
            { "tstop",             LUA_TNUMBER,   NULL,                    TRUE, &tstop },
            { "dt",                LUA_TNUMBER,   NULL,                    TRUE, &dt    },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 4:
            pot = checkPotential(luaSt, 1);
            ic = checkInitialConditions(luaSt, 2);
            tstop = luaL_checknumber(luaSt, 3);
            dt = luaL_checknumber(luaSt, 4);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 4 arguments");
    }

    icNew = reverseOrbit(pot, *ic, tstop, dt);
    pushInitialConditions(luaSt, &icNew);

    return 1;
}

static void registerReverseOrbit(lua_State* luaSt)
{
    lua_pushcfunction(luaSt, luaReverseOrbit);
    lua_setglobal(luaSt, "reverseOrbit");
}

void registerUtilityFunctions(lua_State* luaSt)
{
    registerReverseOrbit(luaSt);
    registerPlummerTimestepIntegral(luaSt);
}

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

    registerNBodyLua(luaSt);

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


