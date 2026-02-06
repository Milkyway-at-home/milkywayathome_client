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
#include "nbody_lua_bfe.h"
#include "nbody_check_params.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"

BFE* checkBFE(lua_State* luaSt, int idx)
{
    return (BFE*) mw_checknamedudata(luaSt, idx, BFE_TYPE);
};

int pushBFE(lua_state* luaSt, const BFE* p)
{
    return pushType(luaSt, BFE_TYPE, sizeof(BFE), (void*) p);
    
};

static const MWEnumAssociation bfeOptions[] =
{
    { "exp",    EXPBFE },
    { "none",   NoBFE  },
    END_MW_ENUM_ASSOCIATION
};

static int createBFE(lua_state* luaSt, const MWNamedArg* argTable, BFE* b)
{
    oneTableArgument(luaSt, argTable);
    if (checkBFEConstants(b))
        luaL_error(luaSt, "Invalid BFE encountered.");
    
    pushBFE(luaSt, b);
    return 1;
};

static int createEXP_BFE(lua_state* luaSt)
{
    static BFE b = EMPTY_BFE;
    static const MWNamedArg argTable[] =
        {
            { "type", LUA_TNUMBER, NULL, TRUE, &b.type, 1},
            { "exp_bfe", LUA_TNUMBER, NULL, TRUE, &b.exp_bfe, 1},
            END_MW_NAMED_ARG
            // not necessarily sure about the notation in line 64
            // pointer to a pointer?
            // will be sorted out later
        };
    
    b.type = EXPBFE;
    return createBFE(luaSt, argTable, &b);
};
    
static int createNo_BFE(lua_state* luaSt)
{
    static BFE b = EMPTY_BFE;
    static const MWNamedArg argTable[] =
        {
            { "type", LUA_TNUMBER, NULL, TRUE, &b.type, 1},
            END_MW_NAMED_ARG
        };
    
    b.type = NoBFE;
    return createBFE(luaSt, argTable, &b);
};

int getBFE_T(lua_state* luaSt, void* v)
{
    return pushEnum(luaSt, bfeOptions, *(int*) v);
};

static int toStringBFE(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showBFE, (LuaTypeCheckFunc) checkBFE);
};

static int eqBFE(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalBFE(checkBFE(luaSt, 1),
                                    checkBFE(luaSt, 2)));
    return 1;
};

int getBFE(lua_State* luaSt, void* v)
{
    pushBFE(luaSt, (BFE*) v);
    return 1;
};

int setSpherical(lua_State* luaSt, void* v)
{
    *(BFE*) v = *checkBFE(luaSt, 2);
    // believe that number corresponds to number of types?
}

static const luaL_reg metaMethodsBFE[] =
{
    { "__tostring", toStringSpherical },
    { "__eq",       eqSpherical       },
    { NULL, NULL }
};

static const luaL_reg methodsBFE[] =
{
    { "EXP_BFE", createEXP_BFE},
    { "none", createNo_BFE},
    { NULL, NULL }
};

static const Xet_reg_pre gettersBFE[] =
{
    { "type", getBFE_T, offsetof(BFE, type) },
    { "exp_bfe", getNumber, offsetof(BFE, exp_bfe)},
    { NULL, NULL, 0}
};

static const Xet_reg_pre settersBFE[] =
{
    { NULL, NULL, 0}
};

int registerBFE(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          BFE_TYPE,
                          gettersBFE,
                          settersBFE,
                          metaMethodsBFE,
                          methodsBFE);
}

int registerBFEKinds(lua_State* luaSt)
{
    int table;
    
    lua_newtable(luaSt);
    table =  lua_gettop(luaSt);
    
    //setModelTableItem(luaSt, table,
}
