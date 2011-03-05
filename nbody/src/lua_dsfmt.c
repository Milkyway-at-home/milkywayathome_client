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

#include "lua_type_marshal.h"
#include "lua_dsfmt.h"
#include "lua_vector.h"

#include "milkyway_util.h"

dsfmt_t* checkDSFMT(lua_State* luaSt, int idx)
{
    return (dsfmt_t*) mw_checknamedudata(luaSt, idx, DSFMT_TYPE);
}

int pushDSFMT(lua_State* luaSt, const dsfmt_t* d)
{
    dsfmt_t* ld;

    ld = (dsfmt_t*) lua_newuserdata(luaSt, sizeof(dsfmt_t));
    if (!ld)
    {
        warn("Creating DSFMT userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, DSFMT_TYPE);
    lua_setmetatable(luaSt, -2);

    *ld = *d;

    return 0;
}

/* TODO: dsfmt_init_by_array */
static int createDSFMT(lua_State* luaSt)
{
    dsfmt_t state;
    uint32_t seed;

    if (lua_gettop(luaSt) > 1)
        return luaL_argerror(luaSt, 1, "Expected 0 or 1 argument");

    /* Omitted or nil argument, seed from clock */
    seed = (uint32_t) luaL_optinteger(luaSt, 1, (uint32_t) time(NULL));

    dsfmt_init_gen_rand(&state, seed);
    pushDSFMT(luaSt, &state);
    return 1;
}

static int dsfmtGenrandOpenOpen(lua_State* luaSt)
{
    lua_pushnumber(luaSt, (lua_Number) dsfmt_genrand_open_open(checkDSFMT(luaSt, 1)));
    return 1;
}

static int dsfmtGenrandClose1Open2(lua_State* luaSt)
{
    lua_pushnumber(luaSt, (lua_Number) dsfmt_genrand_close1_open2(checkDSFMT(luaSt, 1)));
    return 1;
}

static int dsfmtGenrandCloseOpen(lua_State* luaSt)
{
    lua_pushnumber(luaSt, (lua_Number) dsfmt_genrand_close_open(checkDSFMT(luaSt, 1)));
    return 1;
}

static int dsfmtGenrandOpenClose(lua_State* luaSt)
{
    lua_pushnumber(luaSt, (lua_Number) dsfmt_genrand_open_close(checkDSFMT(luaSt, 1)));
    return 1;
}

static int dsfmtRandomVector(lua_State* luaSt)
{
    pushVector(luaSt, mwRandomVector(checkDSFMT(luaSt, 1)));
    return 1;
}

static int dsfmtRandomRange(lua_State* luaSt)
{
    dsfmt_t* d;
    int nArgs;
    double randVal, low, high;

    nArgs = lua_gettop(luaSt);
    d = checkDSFMT(luaSt, 1);

    switch (nArgs - 1)  /* dsfmt_t counts as 1 argument */
    {
        case 0:
            randVal = dsfmt_genrand_open_open(d);
            break;

        case 2:    /* Random number in range */
            low = luaL_checknumber(luaSt, 2);
            high = luaL_checknumber(luaSt, 3);
            randVal = mwXrandom(d, low, high);
            break;
        default:
            return luaL_argerror(luaSt, 1, "Expected 0 or 2 arguments");
    }

    lua_pushnumber(luaSt, (lua_Number) randVal);
    return 1;
}

static int toStringDSFMT(lua_State* luaSt)
{
    lua_pushstring(luaSt, "DSFMT");
    return 1;
}


static const luaL_reg metaMethodsDSFMT[] =
{
    { "__tostring", toStringDSFMT },
    { NULL, NULL }
};

static const luaL_reg methodsDSFMT[] =
{
    { "create",             createDSFMT             },

    /* Lower level bindings */
    { "genrandOpenOpen",    dsfmtGenrandOpenOpen    },
    { "genrandClose1Open2", dsfmtGenrandClose1Open2 },
    { "genrandCloseOpen",   dsfmtGenrandCloseOpen   },
    { "genrandOpenClose",   dsfmtGenrandOpenClose   },

    /* Nicer bindings */
    { "random",             dsfmtRandomRange  },
    { "randomVector",       dsfmtRandomVector },
    { NULL, NULL }
};

static const Xet_reg_pre gettersDSFMT[] =
{
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersDSFMT[] =
{
    { NULL, NULL, 0 }
};

int registerDSFMT(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          DSFMT_TYPE,
                          gettersDSFMT,
                          settersDSFMT,
                          metaMethodsDSFMT,
                          methodsDSFMT);
}

