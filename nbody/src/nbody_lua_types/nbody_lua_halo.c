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
#include "nbody_lua_halo.h"
#include "nbody_check_params.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"

Halo* checkHalo(lua_State* luaSt, int idx)
{
    return (Halo*) mw_checknamedudata(luaSt, idx, HALO_TYPE);
}

int pushHalo(lua_State* luaSt, const Halo* p)
{
    return pushType(luaSt, HALO_TYPE, sizeof(Halo), (void*) p);
}

static const MWEnumAssociation haloOptions[] =
{
    { "logarithmic", LogarithmicHalo },
    { "nfw",         NFWHalo,        },
    { "triaxial",    TriaxialHalo,   },
    END_MW_ENUM_ASSOCIATION
};

static int createHalo(lua_State* luaSt, const MWNamedArg* argTable, Halo* h)
{
    oneTableArgument(luaSt, argTable);
    if (checkHaloConstants(h))
        luaL_error(luaSt, "Invalid halo encountered");

    pushHalo(luaSt, h);
    return 1;
}

static int createLogarithmicHalo(lua_State* luaSt)
{
    static Halo h = EMPTY_HALO;
    static const MWNamedArg argTable[] =
        {
            { "vhalo",       LUA_TNUMBER, NULL, TRUE, &h.vhalo       },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            { "flattenZ",    LUA_TNUMBER, NULL, TRUE, &h.flattenZ    },
            END_MW_NAMED_ARG
        };

    h.type = LogarithmicHalo;
    return createHalo(luaSt, argTable, &h);
}

static int createTriaxialHalo(lua_State* luaSt)
{
    static Halo h = EMPTY_HALO;
    static const MWNamedArg argTable[] =
        {
            { "vhalo",       LUA_TNUMBER, NULL, TRUE, &h.vhalo       },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            { "flattenX",    LUA_TNUMBER, NULL, TRUE, &h.flattenX    },
            { "flattenY",    LUA_TNUMBER, NULL, TRUE, &h.flattenY    },
            { "flattenZ",    LUA_TNUMBER, NULL, TRUE, &h.flattenZ    },
            { "triaxAngle",  LUA_TNUMBER, NULL, TRUE, &h.triaxAngle  },
            END_MW_NAMED_ARG
        };

    h.type = TriaxialHalo;
    return createHalo(luaSt, argTable, &h);
}

static int createNFWHalo(lua_State* luaSt)
{
    static Halo h = EMPTY_HALO;
    static const MWNamedArg argTable[] =
        {
            { "vhalo",       LUA_TNUMBER, NULL, TRUE, &h.vhalo       },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            END_MW_NAMED_ARG
        };

    h.type = NFWHalo;
    return createHalo(luaSt, argTable, &h);
}

static int toStringHalo(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showHalo, (LuaTypeCheckFunc) checkHalo);
}

static int eqHalo(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalHalo(checkHalo(luaSt, 1), checkHalo(luaSt, 2)));
    return 1;
}

int getHaloT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, haloOptions, *(int*) v);
}

int getHalo(lua_State* luaSt, void* v)
{
    pushHalo(luaSt, (Halo*) v);
    return 1;
}

int setHalo(lua_State* luaSt, void* v)
{
    *(Halo*) v = *checkHalo(luaSt, 3);
    return 0;
}

static const luaL_reg metaMethodsHalo[] =
{
    { "__tostring", toStringHalo },
    { "__eq",       eqHalo       },
    { NULL, NULL }
};

static const luaL_reg methodsHalo[] =
{
    { "logarithmic", createLogarithmicHalo },
    { "nfw",         createNFWHalo },
    { "triaxial",    createTriaxialHalo },
    { NULL, NULL }
};

/* TODO Error when writing to fields a halo type doesn't have */
static const Xet_reg_pre gettersHalo[] =
{
    { "type",        getHaloT,  offsetof(Halo, type) },
    { "vhalo",       getNumber, offsetof(Halo, vhalo) },
    { "scaleLength", getNumber, offsetof(Halo, scaleLength) },
    { "flattenX",    getNumber, offsetof(Halo, flattenX) },
    { "flattenY",    getNumber, offsetof(Halo, flattenY) },
    { "flattenZ",    getNumber, offsetof(Halo, flattenZ) },
    { "triaxAngle",  getNumber, offsetof(Halo, triaxAngle) },
    { "c1",          getNumber, offsetof(Halo, c1) },
    { "c2",          getNumber, offsetof(Halo, c2) },
    { "c3",          getNumber, offsetof(Halo, c3) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersHalo[] =
{
    { "vhalo",       setNumber, offsetof(Halo, vhalo) },
    { "scaleLength", setNumber, offsetof(Halo, scaleLength) },
    { "flattenX",    setNumber, offsetof(Halo, flattenX) },
    { "flattenY",    setNumber, offsetof(Halo, flattenY) },
    { "flattenZ",    setNumber, offsetof(Halo, flattenZ) },
    { "triaxAngle",  setNumber, offsetof(Halo, triaxAngle) },
    { "c1",          setNumber, offsetof(Halo, c1) },
    { "c2",          setNumber, offsetof(Halo, c2) },
    { "c3",          setNumber, offsetof(Halo, c3) },
    { NULL, NULL, 0 }
};

int registerHalo(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          HALO_TYPE,
                          gettersHalo,
                          settersHalo,
                          metaMethodsHalo,
                          methodsHalo);
}

/* Add a table with available halo models */
int registerHaloKinds(lua_State* luaSt)
{
    int table;

    lua_newtable(luaSt);
    table = lua_gettop(luaSt);

    setModelTableItem(luaSt, table, createLogarithmicHalo, "logarithmic");
    setModelTableItem(luaSt, table, createNFWHalo, "nfw");
    setModelTableItem(luaSt, table, createTriaxialHalo, "triaxial");

    /* Getting the number of keys in a table is a pain */
    lua_pushnumber(luaSt, 3);
    lua_setfield(luaSt, table, "_count");

    lua_setglobal(luaSt, "haloModels");

    return 0;
}

