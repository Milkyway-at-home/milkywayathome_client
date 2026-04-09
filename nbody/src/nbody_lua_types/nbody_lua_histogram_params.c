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
#include <lauxlib.h>

#include "nbody_types.h"
#include "nbody_show.h"
#include "nbody_lua_histogram_params.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"
#include "nbody_defaults.h"

HistogramParams* checkHistogramParams(lua_State* luaSt, int idx)
{
    return (HistogramParams*) mw_checknamedudata(luaSt, idx, HISTOGRAM_PARAMS_TYPE);
}

int pushHistogramParams(lua_State* luaSt, const HistogramParams* p)
{
    return pushType(luaSt, HISTOGRAM_PARAMS_TYPE, sizeof(HistogramParams), (void*) p);
}

HistogramParams* toHistogramParams(lua_State* luaSt, int idx)
{
    return (HistogramParams*) mw_tonamedudata(luaSt, idx, HISTOGRAM_PARAMS_TYPE);
}

HistogramParams* expectHistogramParams(lua_State* luaSt, int idx)
{
    return (HistogramParams*) expectType(luaSt, idx, HISTOGRAM_PARAMS_TYPE);
}

static int createHistogramParams(lua_State* luaSt)
{
    int nArgs;
    static HistogramParams hp = EMPTY_HISTOGRAM_PARAMS;

    static MWNamedArg argTable[] =
        {
            { "phi",         LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.phi,                1 },
            { "theta",       LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.theta,              1 },
            { "psi",         LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.psi,                1 },
            { "lambdaStart", LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.lambdaStart,        1 },
            { "lambdaEnd",   LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.lambdaEnd,          1 },
            { "lambdaBins",  LUA_TNUMBER,    UINT_TYPE,     TRUE,  &hp.lambdaBins,         1 },
            { "betaStart",   LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.betaStart,          1 },
            { "betaEnd",     LUA_TNUMBER,    REAL_TYPE,     TRUE,  &hp.betaEnd,            1 },
            { "betaBins",    LUA_TNUMBER,    UINT_TYPE,     TRUE,  &hp.betaBins,           1 },
            { "L",           LUA_TTABLE,     REAL_TYPE,     FALSE, &hp.L,                  3 }, // These two could probably be input as vectors, but when i tried it didnt work
            { "LErr",        LUA_TTABLE,     REAL_TYPE,     FALSE, &hp.LErr,               3 }, // Using tables should be fine
            { "nRange",      LUA_TNUMBER,    UINT_TYPE,     FALSE, &hp.nRange,             1 },
            { "EMDRange",    LUA_TTABLE,     REAL_TYPE,     FALSE, &hp.EMDRange,           0 }, // This is an array, length set below. Used 0 so ignored if nRange isnt set
            END_MW_NAMED_ARG
        };

    /* Read nRange first to set EMDRange length*/
    unsigned int temp = 0;
    if (lua_istable(luaSt, 1))
    {
        lua_getfield(luaSt, 1, "nRange");
        if(lua_isnil(luaSt, -1))
        {
            lua_pop(luaSt, 1);
        }
        else if(!lua_isnumber(luaSt, -1))
        {
            argTable[12].arrayLen = temp;
        }
        else
        {
            temp = (unsigned int)lua_tonumber(luaSt, -1);
            lua_pop(luaSt, 1);
            argTable[12].arrayLen = temp;
        }
    }

    nArgs = lua_gettop(luaSt);
    if (nArgs == 0)
    {
        pushHistogramParams(luaSt, &defaultHistogramParams);
    }
    else if (nArgs == 1)
    {
        handleNamedArgumentTable(luaSt, argTable, 1);
        pushHistogramParams(luaSt, &hp);
    }
    else
    {
        return luaL_argerror(luaSt, 1, "Expected argument table");
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfloat-equal"
    if(hp.nRange >= 2)
    {
        printf("NOTE: Using EMD Ranges given in lua. Histogram input will be ignored.\n");
        if(hp.EMDRange[hp.nRange - 1] == 0)
        {
            printf("WARNING: EMDRange ends in zero. Is nRange too large?\n");
        }
    }
    #pragma GCC diagnostic pop
    if(hp.L.x > 0.00001 || hp.L.y > 0.00001 || hp.L.z > 0.00001 || hp.L.x < -0.00001 || hp.L.y < -0.00001 || hp.L.z < -0.00001)
    {
        printf("NOTE: Using Momentum info given in lua. Histogram input will be ignored.\n");
    }
    if(hp.nRange % 2 != 0)
    {
        printf("WARNING: nRange is odd, but should be even. Will use nRange - 1.\n");
    }

    return 1;
}

static int toStringHistogramParams(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showHistogramParams, (LuaTypeCheckFunc) checkHistogramParams);
}

static int eqHistogramParams(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalHistogramParams(checkHistogramParams(luaSt, 1), checkHistogramParams(luaSt, 2)));
    return 1;
}

int getHistogramParams(lua_State* luaSt, void* v)
{
    pushHistogramParams(luaSt, (HistogramParams*) v);
    return 1;
}

int setHistogramParams(lua_State* luaSt, void* v)
{
    *(HistogramParams*) v = *checkHistogramParams(luaSt, 3);
    return 0;
}

static const luaL_reg metaMethodsHistogramParams[] =
{
    { "__tostring", toStringHistogramParams },
    { "__eq",       eqHistogramParams       },
    { NULL, NULL }
};

static const luaL_reg methodsHistogramParams[] =
{
    { "create", createHistogramParams },
    { NULL, NULL }
};

int getEMDRange(lua_State* luaSt, void* v) 
{
    HistogramParams* hp = (HistogramParams*)((char*)v - offsetof(HistogramParams, EMDRange));
    return getRealArray(luaSt, hp->EMDRange, hp->nRange);
}

static const Xet_reg_pre gettersHistogramParams[] =
{
    { "phi" ,        getNumber,    offsetof(HistogramParams, phi)          },
    { "theta",       getNumber,    offsetof(HistogramParams, theta)        },
    { "psi",         getNumber,    offsetof(HistogramParams, psi)          },
    { "lambdaStart", getNumber,    offsetof(HistogramParams, lambdaStart)  },
    { "lambdaEnd",   getNumber,    offsetof(HistogramParams, lambdaEnd)    },
    { "lambdaBins",  getUInt,      offsetof(HistogramParams, lambdaBins)   },
    { "betaStart",   getNumber,    offsetof(HistogramParams, betaStart)    },
    { "betaEnd",     getNumber,    offsetof(HistogramParams, betaEnd)      },
    { "betaBins",    getUInt,      offsetof(HistogramParams, betaBins)     },
    { "L",           getVector,    offsetof(HistogramParams, L)            }, 
    { "LErr",        getVector,    offsetof(HistogramParams, LErr)         },
    { "nRange",      getUInt,      offsetof(HistogramParams, nRange)       },
    { "EMDRange",    getEMDRange,  offsetof(HistogramParams, EMDRange)     },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersHistogramParams[] =
{
    { "phi" ,        setNumber,    offsetof(HistogramParams, phi)         },
    { "theta",       setNumber,    offsetof(HistogramParams, theta)       },
    { "psi",         setNumber,    offsetof(HistogramParams, psi)         },
    { "lambdaStart", setNumber,    offsetof(HistogramParams, lambdaStart) },
    { "lambdaEnd",   setNumber,    offsetof(HistogramParams, lambdaEnd)   },
    { "lambdaBins",  setUInt,      offsetof(HistogramParams, lambdaBins)  },
    { "betaStart",   setNumber,    offsetof(HistogramParams, betaStart)   },
    { "betaEnd",     setNumber,    offsetof(HistogramParams, betaEnd)     },
    { "betaBins",    setUInt,      offsetof(HistogramParams, betaBins)    },
    { "L",           setVector,    offsetof(HistogramParams, L)           },
    { "LErr",        setVector,    offsetof(HistogramParams, LErr)        },
    { "nRange",      setUInt,      offsetof(HistogramParams, nRange)      },
    { "EMDRange",    setRealArray, offsetof(HistogramParams, EMDRange)    },
    { NULL, NULL, 0 }
};

int registerHistogramParams(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          HISTOGRAM_PARAMS_TYPE,
                          gettersHistogramParams,
                          settersHistogramParams,
                          metaMethodsHistogramParams,
                          methodsHistogramParams);
}

