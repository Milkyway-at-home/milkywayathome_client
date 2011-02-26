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
#include "show.h"

#include "lua_type_marshal.h"
#include "lua_vector.h"
#include "lua_initial_conditions.h"
#include "lua_nbodyctx.h"

#include "milkyway_util.h"
#include "nbody_util.h"


InitialConditions* checkInitialConditions(lua_State* luaSt, int index)
{
    return (InitialConditions*) mw_checknamedudata(luaSt, index, INITIAL_CONDITIONS_TYPE);
}

int pushInitialConditions(lua_State* luaSt, const InitialConditions* ic)
{
    InitialConditions* lic;

    lic = (InitialConditions*)lua_newuserdata(luaSt, sizeof(InitialConditions));
    if (!lic)
    {
        warn("Creating InitialConditions userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, INITIAL_CONDITIONS_TYPE);
    lua_setmetatable(luaSt, -2);

    *lic = *ic;

    return 0;
}


/* Arguments: InitialConditions.create(NBodyCtx, position, velocity, [useGalacticCoordinates, useRadians])
   or named arguments, e.g. InitialConditions.create{ context = ctx,
                                                      position = x,
                                                      velocity = v,
                                                      useGalacticCoordinates = false,
                                                      useRadians = false)
 */
static int createInitialConditions(lua_State* luaSt)
{
    int nArgs;
    const NBodyCtx* ctx = NULL;
    mwbool useGalacticCoordinates = FALSE;
    mwbool useRadians = FALSE;
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;

    const MWNamedArg argTable[] =
        {
            { "context",                LUA_TUSERDATA, NBODY_CTX, TRUE,  &ctx                    },
            { "useGalacticCoordinates", LUA_TBOOLEAN,  NULL,      FALSE, &useGalacticCoordinates },
            { "useRadians",             LUA_TBOOLEAN,  NULL,      FALSE, &useRadians             },
            END_MW_NAMED_ARG
        };

    warn("Creating initial conditions\n");

    nArgs = lua_gettop(luaSt);
    switch (nArgs)
    {
        case 1:  /* Table of named arguments */
            /* TODO */
            handleNamedArgumentTable(luaSt, argTable, 1);

            //mw_panic("Implement me!\n");
            break;

        case 3:  /* (Context, Position, velocity) Use defaults for rest */
            ctx = checkNBodyCtx(luaSt, 1);
            ic.position = *checkVector(luaSt, 2);
            ic.velocity = *checkVector(luaSt, 3);
            break;

        case 5: /* Everything */
            ctx = checkNBodyCtx(luaSt, 1);
            ic.position = *checkVector(luaSt, 2);
            ic.velocity = *checkVector(luaSt, 3);
            useGalacticCoordinates = mw_lua_checkboolean(luaSt, 4);
            useRadians = mw_lua_checkboolean(luaSt, 5);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1, 3, or 5 arguments");

    }


    warn("SLURP\n");
    printNBodyCtx(ctx);

    /* We aren't given galactic coordinates, so convert them */
    if (!useGalacticCoordinates)
    {
        ic.position = useRadians ? lbrToCartesian_rad(ctx, ic.position) : lbrToCartesian(ctx, ic.position);
    }

    pushInitialConditions(luaSt, &ic);
    return 1;
}

static int toStringInitialConditions(lua_State* luaSt)
{
    InitialConditions* ic;
    char* str;

    ic = checkInitialConditions(luaSt, 1);
    str = showInitialConditions(ic);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static const luaL_reg metaMethodsInitialConditions[] =
{
    { "__tostring", toStringInitialConditions },
    { NULL, NULL }
};

static const luaL_reg methodsInitialConditions[] =
{
    { "create",   createInitialConditions },
    { NULL, NULL }
};

/* TODO Error when writing to fields a halo type doesn't have */
static const Xet_reg_pre gettersInitialConditions[] =
{
    { "position",     getVector, offsetof(InitialConditions, position)     },
    { "velocity",     getVector, offsetof(InitialConditions, velocity)     },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersInitialConditions[] =
{
    { "position",     setVector, offsetof(InitialConditions, position)     },
    { "velocity",     setVector, offsetof(InitialConditions, velocity)     },
    { NULL, NULL, 0 }
};

int registerInitialConditions(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          INITIAL_CONDITIONS_TYPE,
                          gettersInitialConditions,
                          settersInitialConditions,
                          metaMethodsInitialConditions,
                          methodsInitialConditions);
}

