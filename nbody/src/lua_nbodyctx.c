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
#include "lua_type_marshal.h"
#include "lua_nbodyctx.h"
#include "show.h"

#include "milkyway_util.h"

static const NBodyCtx defaultNBodyCtx =
{
    /* Grr lack of C99 named struct initializers in MSVC */
    /* .pot             */  EMPTY_POTENTIAL,
    /* .nbody           */  0,
    /* .timestep        */  0.0,
    /* .timeEvolve      */  0.0,

    /* .theta           */  0.0,
    /* .eps2            */  0.0,

    /* .treeRSize       */  4.0,
    /* .sunGCDist       */  8.0,
    /* .criterion       */  NewCriterion,
    /* .useQuad         */  TRUE,
    /* .allowIncest     */  FALSE,

    /* .outputCartesian */  FALSE,
    /* .outputBodies    */  FALSE,
    /* .outputHistogram */  FALSE,

    /* .histogramParams */  EMPTY_HISTOGRAM_PARAMS,

    /* .freqout         */  0,
    /* .outfile         */  NULL
};


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
    warn("Set criterion\n");
    *(criterion_t*) v = checkEnum(luaSt, criterionOptions, -1);
    return 0;
}

NBodyCtx* checkNBodyCtx(lua_State* luaSt, int idx)
{
    return (NBodyCtx*) mw_checknamedudata(luaSt, idx, NBODY_CTX);
}

int pushNBodyCtx(lua_State* luaSt, const NBodyCtx* ctx)
{
    NBodyCtx* lctx;

    lctx = (NBodyCtx*)lua_newuserdata(luaSt, sizeof(NBodyCtx));

    luaL_getmetatable(luaSt, NBODY_CTX);
    lua_setmetatable(luaSt, -2);

    *lctx = *ctx;

    return 0;
}

static int createNBodyCtx(lua_State* luaSt)
{
    static NBodyCtx ctx;
    static real freqOutf = 0.0;
    static const char* criterionName = NULL;

    static const MWNamedArg argTable[] =
        {
            { "timestep",    LUA_TNUMBER,  NULL, TRUE,  &ctx.timestep    },
            { "timeEvolve",  LUA_TNUMBER,  NULL, TRUE,  &ctx.timeEvolve  },
            { "freqOut",     LUA_TNUMBER,  NULL, FALSE, &freqOutf        },
            { "theta",       LUA_TNUMBER,  NULL, TRUE,  &ctx.theta       },
            { "eps2",        LUA_TNUMBER,  NULL, TRUE,  &ctx.eps2        },
            { "treeRSize",   LUA_TNUMBER,  NULL, FALSE, &ctx.treeRSize   },
            { "sunGCDist",   LUA_TNUMBER,  NULL, FALSE, &ctx.sunGCDist   },
            { "criterion",   LUA_TSTRING,  NULL, FALSE, &criterionName   },
            { "useQuad",     LUA_TBOOLEAN, NULL, FALSE, &ctx.useQuad     },
            { "allowIncest", LUA_TBOOLEAN, NULL, FALSE, &ctx.allowIncest },
            END_MW_NAMED_ARG
        };

    criterionName = NULL;
    freqOutf = 0.0;
    ctx = defaultNBodyCtx;

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected named argument table");

    handleNamedArgumentTable(luaSt, argTable, 1);

    ctx.freqOut = (unsigned int) freqOutf;

    /* FIXME: Hacky handling of enum. Will result in not good error
     * messages as well as not fitting in. */
    if (criterionName) /* Not required */
    {
        lua_pushstring(luaSt, criterionName);
        ctx.criterion = checkEnum(luaSt, criterionOptions, lua_gettop(luaSt));
        lua_pop(luaSt, 1);
    }

    pushNBodyCtx(luaSt, &ctx);
    return 1;
}

static int toStringNBodyCtx(lua_State* luaSt)
{
    NBodyCtx* ctx;
    char* str;

    ctx = checkNBodyCtx(luaSt, 1);
    str = showNBodyCtx(ctx);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static const luaL_reg metaMethodsNBodyCtx[] =
{
    { "__tostring", toStringNBodyCtx },
    { NULL, NULL }
};

static const luaL_reg methodsNBodyCtx[] =
{
    { "create", createNBodyCtx },
    { NULL, NULL }
};

static const Xet_reg_pre gettersNBodyCtx[] =
{
    { "timestep",        getNumber,     offsetof(NBodyCtx, timestep)        },
    { "timeEvolve",      getNumber,     offsetof(NBodyCtx, timeEvolve)      },
    { "freqOut",         getNumber,     offsetof(NBodyCtx, freqOut)         },
    { "theta",           getNumber,     offsetof(NBodyCtx, theta)           },
    { "eps2",            getNumber,     offsetof(NBodyCtx, eps2)            },
    { "treeRSize",       getNumber,     offsetof(NBodyCtx, treeRSize)       },
    { "sunGCDist",       getNumber,     offsetof(NBodyCtx, sunGCDist)       },
    { "criterion",       getCriterionT, offsetof(NBodyCtx, criterion)       },
    { "useQuad",         getBool,       offsetof(NBodyCtx, useQuad)         },
    { "allowIncest",     getBool,       offsetof(NBodyCtx, allowIncest)     },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersNBodyCtx[] =
{
    { "timestep",        setNumber,     offsetof(NBodyCtx, timestep)        },
    { "timeEvolve",      setNumber,     offsetof(NBodyCtx, timeEvolve)      },
    { "freqOut",         setNumber,     offsetof(NBodyCtx, freqOut)         },
    { "theta",           setNumber,     offsetof(NBodyCtx, theta)           },
    { "eps2",            setNumber,     offsetof(NBodyCtx, eps2)            },
    { "treeRSize",       setNumber,     offsetof(NBodyCtx, treeRSize)       },
    { "sunGCDist",       setNumber,     offsetof(NBodyCtx, sunGCDist)       },
    { "criterion",       setCriterionT, offsetof(NBodyCtx, criterion)       },
    { "useQuad",         setBool,       offsetof(NBodyCtx, useQuad)         },
    { "allowIncest",     setBool,       offsetof(NBodyCtx, allowIncest)     },
    { NULL, NULL, 0 }
};

int registerNBodyCtx(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          NBODY_CTX,
                          gettersNBodyCtx,
                          settersNBodyCtx,
                          metaMethodsNBodyCtx,
                          methodsNBodyCtx);
}

