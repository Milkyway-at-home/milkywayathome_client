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

    static const MWNamedArg argTable[] =
        {
            { "phi",      LUA_TNUMBER, NULL, TRUE, &hp.phi      },
            { "theta",    LUA_TNUMBER, NULL, TRUE, &hp.theta    },
            { "startRaw", LUA_TNUMBER, NULL, TRUE, &hp.startRaw },
            { "endRaw",   LUA_TNUMBER, NULL, TRUE, &hp.endRaw   },
            { "binSize",  LUA_TNUMBER, NULL, TRUE, &hp.binSize  },
            { "center",   LUA_TNUMBER, NULL, TRUE, &hp.center   },
            END_MW_NAMED_ARG
        };

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

static const Xet_reg_pre gettersHistogramParams[] =
{
    { "phi" ,     getNumber, offsetof(HistogramParams, phi)      },
    { "theta",    getNumber, offsetof(HistogramParams, theta)    },
    { "startRaw", getNumber, offsetof(HistogramParams, startRaw) },
    { "endRaw",   getNumber, offsetof(HistogramParams, endRaw)   },
    { "binSize",  getNumber, offsetof(HistogramParams, binSize)  },
    { "center",   getNumber, offsetof(HistogramParams, center)   },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersHistogramParams[] =
{
    { "phi" ,     setNumber, offsetof(HistogramParams, phi)      },
    { "theta",    setNumber, offsetof(HistogramParams, theta)    },
    { "startRaw", setNumber, offsetof(HistogramParams, startRaw) },
    { "endRaw",   setNumber, offsetof(HistogramParams, endRaw)   },
    { "binSize",  setNumber, offsetof(HistogramParams, binSize)  },
    { "center",   setNumber, offsetof(HistogramParams, center)   },
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

