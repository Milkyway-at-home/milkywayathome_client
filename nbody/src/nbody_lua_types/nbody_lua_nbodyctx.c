/*
 * Copyright (c) 2011 Matthew Arsenault
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "nbody_show.h"
#include "nbody_lua_nbodyctx.h"
#include "nbody_lua_potential.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"
#include "nbody_defaults.h"
#include "nbody_lua_misc.h"
#include "nbody_util.h"
#include "nbody_lua_util.h"


static const MWEnumAssociation criterionOptions[] =
{
    { "NewCriterion", NewCriterion },
    { "Exact",        Exact        },
    { "BH86",         BH86         },
    { "SW93",         SW93         },
    END_MW_ENUM_ASSOCIATION
};

static int getCriterionT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, criterionOptions, *(criterion_t*) v);
}

static int setCriterionT(lua_State* luaSt, void* v)
{
    *(criterion_t*) v = checkEnum(luaSt, criterionOptions, 3);
    return 0;
}

NBodyCtx* checkNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) mw_checknamedudata(luaSt, idx, NBODYCTX_TYPE);
}

NBodyCtx* toNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) mw_tonamedudata(luaSt, idx, NBODYCTX_TYPE);
}

NBodyCtx* expectNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) expectType(luaSt, idx, NBODYCTX_TYPE);
}

int pushNBodyCtx(lua_State* luaSt, const NBodyCtx* p)
{
    return pushType(luaSt, NBODYCTX_TYPE, sizeof(NBodyCtx), (void*) p);
}

criterion_t readCriterion(lua_State* luaSt, const char* name)
{
    return (criterion_t) readEnum(luaSt, criterionOptions, name);
}

static int createNBodyCtx(lua_State* luaSt)
{
    static NBodyCtx ctx;
    static const char* criterionName = NULL;
    real nStepf = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "timestep",    LUA_TNUMBER,  NULL, TRUE,  &ctx.timestep    },
            { "timeEvolve",  LUA_TNUMBER,  NULL, TRUE,  &ctx.timeEvolve  },
            { "theta",       LUA_TNUMBER,  NULL, FALSE, &ctx.theta       },
            { "eps2",        LUA_TNUMBER,  NULL, TRUE,  &ctx.eps2        },
            { "treeRSize",   LUA_TNUMBER,  NULL, FALSE, &ctx.treeRSize   },
            { "sunGCDist",   LUA_TNUMBER,  NULL, FALSE, &ctx.sunGCDist   },
            { "criterion",   LUA_TSTRING,  NULL, FALSE, &criterionName   },
            { "useQuad",     LUA_TBOOLEAN, NULL, FALSE, &ctx.useQuad     },
            { "allowIncest", LUA_TBOOLEAN, NULL, FALSE, &ctx.allowIncest },
            { "quietErrors", LUA_TBOOLEAN, NULL, FALSE, &ctx.quietErrors },
            END_MW_NAMED_ARG
        };

    criterionName = NULL;
    ctx = defaultNBodyCtx;

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected named argument table");

    handleNamedArgumentTable(luaSt, argTable, 1);

    /* FIXME: Hacky handling of enum. Will result in not good error
     * messages as well as not fitting in. */
    if (criterionName) /* Not required */
    {
        ctx.criterion = readCriterion(luaSt, criterionName);
    }

    if ((ctx.criterion != Exact) && (ctx.theta < 0.0))
    {
        return luaL_argerror(luaSt, 1, "Theta argument required for criterion != 'Exact'");
    }
    else if (ctx.criterion == Exact)
    {
        /* These don't mean anything here */
        ctx.theta = 0.0;
        ctx.useQuad = FALSE;
    }

    nStepf = mw_ceil(ctx.timeEvolve / ctx.timestep);
    if (nStepf >= (real) UINT_MAX)
    {
        luaL_error(luaSt,
                   "Number of timesteps exceeds UINT_MAX: %f timesteps (%f / %f)\n",
                   nStepf,
                   ctx.timeEvolve, ctx.timestep);
    }

    ctx.nStep = (unsigned int) nStepf;

    {
        int major = 0, minor = 0;

        /* Automatically correct the timestep size so an integer
         * number of timesteps covers the evolution time.
         *
         * Only do this if we require a minimum version of 0.90 to
         * avoid not validating against currently existing workunits
         */

        if (  !nbReadMinVersion(luaSt, &major, &minor)    /* If we fail to read version */
            || (major > 0 || (major == 0 && minor >= 90)) /* Version required >= 0.90 */
            || (major == 0 && minor == 0))                /* Min version not set */
        {
            ctx.timestep = nbCorrectTimestep(ctx.timeEvolve, ctx.timestep);
        }
        else
        {
            mw_printf("Warning: not applying timestep correction for workunit with min version %d.%d\n", major, minor);
        }
    }

    pushNBodyCtx(luaSt, &ctx);
    return 1;
}

static int toStringNBodyCtx(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showNBodyCtx, (LuaTypeCheckFunc) checkNBodyCtx);
}

static int eqNBodyCtx(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalNBodyCtx(checkNBodyCtx(luaSt, 1), checkNBodyCtx(luaSt, 2)));
    return 1;
}

static int addPotential(lua_State* luaSt)
{
    NBodyCtx* ctx;

    if (lua_gettop(luaSt) != 2)
        return luaL_argerror(luaSt, 1, "Expected named 2 arguments");

    ctx = checkNBodyCtx(luaSt, 1);
    ctx->pot = *checkPotential(luaSt, 2);

    return 0;
}

static const luaL_reg metaMethodsNBodyCtx[] =
{
    { "__tostring", toStringNBodyCtx },
    { "__eq",       eqNBodyCtx       },
    { NULL, NULL }
};

static const luaL_reg methodsNBodyCtx[] =
{
    { "create", createNBodyCtx },
    { "addPotential", addPotential },
    { NULL, NULL }
};

static const Xet_reg_pre gettersNBodyCtx[] =
{
    { "timestep",        getNumber,     offsetof(NBodyCtx, timestep)    },
    { "timeEvolve",      getNumber,     offsetof(NBodyCtx, timeEvolve)  },
    { "theta",           getNumber,     offsetof(NBodyCtx, theta)       },
    { "eps2",            getNumber,     offsetof(NBodyCtx, eps2)        },
    { "treeRSize",       getNumber,     offsetof(NBodyCtx, treeRSize)   },
    { "sunGCDist",       getNumber,     offsetof(NBodyCtx, sunGCDist)   },
    { "criterion",       getCriterionT, offsetof(NBodyCtx, criterion)   },
    { "useQuad",         getBool,       offsetof(NBodyCtx, useQuad)     },
    { "allowIncest",     getBool,       offsetof(NBodyCtx, allowIncest) },
    { "quietErrors",     getBool,       offsetof(NBodyCtx, quietErrors) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersNBodyCtx[] =
{
    { "timestep",        setNumber,     offsetof(NBodyCtx, timestep)    },
    { "timeEvolve",      setNumber,     offsetof(NBodyCtx, timeEvolve)  },
    { "theta",           setNumber,     offsetof(NBodyCtx, theta)       },
    { "eps2",            setNumber,     offsetof(NBodyCtx, eps2)        },
    { "treeRSize",       setNumber,     offsetof(NBodyCtx, treeRSize)   },
    { "sunGCDist",       setNumber,     offsetof(NBodyCtx, sunGCDist)   },
    { "criterion",       setCriterionT, offsetof(NBodyCtx, criterion)   },
    { "useQuad",         setBool,       offsetof(NBodyCtx, useQuad)     },
    { "allowIncest",     setBool,       offsetof(NBodyCtx, allowIncest) },
    { "quietErrors",     setBool,       offsetof(NBodyCtx, quietErrors) },
    { NULL, NULL, 0 }
};

int registerNBodyCtx(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          NBODYCTX_TYPE,
                          gettersNBodyCtx,
                          settersNBodyCtx,
                          metaMethodsNBodyCtx,
                          methodsNBodyCtx);
}

