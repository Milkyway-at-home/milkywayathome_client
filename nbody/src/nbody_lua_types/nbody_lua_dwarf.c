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
    static Dwarf h = EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass       },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            END_MW_NAMED_ARG
        };

    h.type = Plummer;
    return createDwarf(luaSt, argTable, &h);
}

static int createNFWDwarf(lua_State* luaSt)
{
    static Dwarf h = EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass       },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            END_MW_NAMED_ARG
        };

    h.type = NFW;
    return createDwarf(luaSt, argTable, &h);
}

static int createGen_HernDwarf(lua_State* luaSt)
{
    static Dwarf h = EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass       },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            END_MW_NAMED_ARG
        };

    h.type = General_Hernquist;
    return createDwarf(luaSt, argTable, &h);
}

static int createEinastoDwarf(lua_State* luaSt)
{
    static Dwarf h = EMPTY_DWARF;
    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &h.mass        },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &h.scaleLength },
            { "n"          , LUA_TNUMBER, NULL, TRUE, &h.n           }, 
            END_MW_NAMED_ARG
        };

    h.type = Einasto;
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

    /* Getting the number of keys in a table is a pain */
    lua_pushnumber(luaSt, 3);
    lua_setfield(luaSt, table, "_count");

    lua_setglobal(luaSt, "dwarfModels");

    return 0;
}

