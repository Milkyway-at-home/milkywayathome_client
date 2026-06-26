/*
Copyright (C) 2011  Matthew Arsenault
Copyright (c) 2016 Siddhartha Shelton
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
#include "nbody_lua_dwarf.h"
#include "nbody_check_params.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"

Dwarf* checkDwarf(lua_State* luaSt, int idx)
{
    return (Dwarf*) mw_checknamedudata(luaSt, idx, DWARF_TYPE);
}

int pushDwarf(lua_State* luaSt, const Dwarf* p)
{
    return pushType(luaSt, DWARF_TYPE, sizeof(Dwarf), (void*) p);
}

static const MWEnumAssociation dwarfOptions[] =
{
    { "plummer",              Plummer              },
    { "nfw",                  NFW,                 },
    { "general_hernquist",    General_Hernquist,   },
    { "einasto",              Einasto,             },
    { "cored",                Cored,               },
    { "king",                 King,                },
    END_MW_ENUM_ASSOCIATION
};

static int createDwarf(lua_State* luaSt, const MWNamedArg* argTable, Dwarf* h)
{
    oneTableArgument(luaSt, argTable);
    pushDwarf(luaSt, h);
    return 1;
}

static int createPlummerDwarf(lua_State* luaSt)
{
    static Dwarf h;
    h = (Dwarf)EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass,        1 },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength, 1 },
            END_MW_NAMED_ARG
        };

    h.type = Plummer;
    oneTableArgument(luaSt, argTable);
    if (h.mass < 0.0)
        luaL_error(luaSt, "Plummer dwarf mass must be non-negative");
    if (h.scaleLength <= 0.0)
        luaL_error(luaSt, "Plummer dwarf scaleLength must be positive");
    return createDwarf(luaSt, argTable, &h);
}

static int createNFWDwarf(lua_State* luaSt)
{
    static Dwarf h;
    h = (Dwarf)EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass,        1 },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength, 1 },
            { "rcut",        LUA_TNUMBER, NULL, FALSE, &h.rcut,       1 },
            END_MW_NAMED_ARG
        };
    
    /* Defaults: rcut = 0.0 (no cutoff) */
    h.type = NFW;
    h.rcut = 0.0;
    oneTableArgument(luaSt, argTable);
    if (h.mass < 0.0)
        luaL_error(luaSt, "NFW dwarf mass must be non-negative");
    if (h.scaleLength <= 0.0)
        luaL_error(luaSt, "NFW dwarf scaleLength must be positive");
    if (h.rcut < 0.0)
        luaL_error(luaSt, "NFW dwarf rcut must be non-negative");
    return createDwarf(luaSt, argTable, &h);
}

static int createGen_HernDwarf(lua_State* luaSt)
{
    static Dwarf h;
    h = (Dwarf)EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass,        1 },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength, 1 },
            END_MW_NAMED_ARG
        };

    h.type = General_Hernquist;
    oneTableArgument(luaSt, argTable);
    if (h.mass < 0.0)
        luaL_error(luaSt, "General Hernquist dwarf mass must be non-negative");
    if (h.scaleLength <= 0.0)
        luaL_error(luaSt, "General Hernquist dwarf scaleLength must be positive");
    return createDwarf(luaSt, argTable, &h);
}

static int createEinastoDwarf(lua_State* luaSt)
{
    static Dwarf h;
    h = (Dwarf)EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass,        1 },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength, 1 },
            { "n"          , LUA_TNUMBER, NULL, TRUE, &h.n,           1 }, 
            END_MW_NAMED_ARG
        };

    h.type = Einasto;
    oneTableArgument(luaSt, argTable);
    if (h.mass < 0.0)
        luaL_error(luaSt, "Einasto dwarf mass must be non-negative");
    if (h.scaleLength <= 0.0)
        luaL_error(luaSt, "Einasto dwarf scaleLength must be positive");
    if (h.n <= 0.0)
        luaL_error(luaSt, "Einasto dwarf n must be positive");
    return createDwarf(luaSt, argTable, &h);
}

static int createCoredDwarf(lua_State* luaSt)
{
    static Dwarf h;
    h = (Dwarf)EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass,        1 },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength, 1 },
			{ "r1", 		 LUA_TNUMBER, NULL, TRUE, &h.r1,          1 },
			{ "rc", 		 LUA_TNUMBER, NULL, TRUE, &h.rc,          1 },
			{ "rcut", 	     LUA_TNUMBER, NULL, FALSE, &h.rcut,       1 },
            END_MW_NAMED_ARG
        };

    h.type = Cored;
	h.r1 = h.scaleLength;
	h.rcut = 0.0;
    oneTableArgument(luaSt, argTable);
    if (h.mass < 0.0)
        luaL_error(luaSt, "Cored dwarf mass must be non-negative");
    if (h.scaleLength <= 0.0)
        luaL_error(luaSt, "Cored dwarf scaleLength must be positive");
    if (h.r1 <= 0.0)
        luaL_error(luaSt, "Cored dwarf r1 must be positive");
    if (h.rc <= 0.0)
        luaL_error(luaSt, "Cored dwarf rc must be positive");
    if (h.rcut < 0.0)
        luaL_error(luaSt, "Cored dwarf rcut must be non-negative");
    if ((mw_abs(h.rcut) > 0.00001) && (h.rcut < h.r1))
        {
        luaL_error(luaSt, "Cored dwarf rcut must be no less than r1");
        }
    return createDwarf(luaSt, argTable, &h);
}

static int createKingDwarf(lua_State* luaSt)
{
    static Dwarf h;
    h = (Dwarf)EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass,        1 },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength, 1 },
			{ "W0", 		 LUA_TNUMBER, NULL, TRUE, &h.W0,          1 },
            END_MW_NAMED_ARG
        };

    h.type = King;
    oneTableArgument(luaSt, argTable);
    if (h.mass < 0.0)
        luaL_error(luaSt, "King model mass must be non-negative");
    if (h.scaleLength <= 0.0)
        luaL_error(luaSt, "King model scaleLength must be positive");
    if (h.W0 <= 0.0)
        luaL_error(luaSt, "King model W0 must be positive");
    return createDwarf(luaSt, argTable, &h);
}


int getDwarfT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, dwarfOptions, *(int*) v);
}

int getDwarf(lua_State* luaSt, void* v)
{
    pushDwarf(luaSt, (Dwarf*) v);
    return 1;
}

int setDwarf(lua_State* luaSt, void* v)
{
    *(Dwarf*) v = *checkDwarf(luaSt, 3);
    return 0;
}

static const luaL_reg metaMethodsDwarf[] =
{
//     { "__tostring", toStringDwarf },
//     { "__eq",       eqDwarf       },
    { NULL, NULL }
};

static const luaL_reg methodsDwarf[] =
{
    { "plummer",              createPlummerDwarf    },
    { "nfw",                  createNFWDwarf        },
    { "general_hernquist",    createGen_HernDwarf   },
    { "einasto",              createEinastoDwarf    },
    { "cored",                createCoredDwarf     },
    { "king",                 createKingDwarf      },
    { NULL, NULL }
};

/* TODO Error when writing to fields a Dwarf type doesn't have */
static const Xet_reg_pre gettersDwarf[] =
{
    { "type",        getDwarfT, offsetof(Dwarf, type) },
    { "mass",        getNumber, offsetof(Dwarf, mass) },
    { "scaleLength", getNumber, offsetof(Dwarf, scaleLength) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersDwarf[] =
{
    { "mass",        setNumber, offsetof(Dwarf, mass) },
    { "scaleLength", setNumber, offsetof(Dwarf, scaleLength) },
    { NULL, NULL, 0 }
};

int registerDwarf(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          DWARF_TYPE,
                          gettersDwarf,
                          settersDwarf,
                          metaMethodsDwarf,
                          methodsDwarf);
}

/* Add a table with available Dwarf models */
int registerDwarfKinds(lua_State* luaSt)
{
    int table;

    lua_newtable(luaSt);
    table = lua_gettop(luaSt);

    setModelTableItem(luaSt, table, createPlummerDwarf, "plummer");
    setModelTableItem(luaSt, table, createNFWDwarf, "nfw");
    setModelTableItem(luaSt, table, createGen_HernDwarf, "general_hernquist");
    setModelTableItem(luaSt, table, createEinastoDwarf, "einasto");
    setModelTableItem(luaSt, table, createCoredDwarf, "cored");
    setModelTableItem(luaSt, table, createKingDwarf, "king");
    
    /* Getting the number of keys in a table is a pain */
    lua_pushnumber(luaSt, 3);
    lua_setfield(luaSt, table, "_count");

    lua_setglobal(luaSt, "dwarfModels");

    return 0;
}

