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
#include "nbody_lua_spherical.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"
#include "nbody_check_params.h"

Spherical* checkSpherical(lua_State* luaSt, int idx)
{
    return (Spherical*) mw_checknamedudata(luaSt, idx, SPHERICAL_TYPE);
}

int pushSpherical(lua_State* luaSt, const Spherical* p)
{
    return pushType(luaSt, SPHERICAL_TYPE, sizeof(Spherical), (void*) p);
}

static const MWEnumAssociation sphericalOptions[] =
{
    { "hernquist", HernquistSpherical },
    { "plummer", PlummerSpherical },
    { "none", NoSpherical },
    END_MW_ENUM_ASSOCIATION
};

static int createSpherical(lua_State* luaSt, const MWNamedArg* argTable, Spherical* s)
{
    oneTableArgument(luaSt, argTable);
    if (checkSphericalConstants(s))
        luaL_error(luaSt, "Invalid bulge encountered");

    pushSpherical(luaSt, s);
    return 1;
}

static int createHernquistSpherical(lua_State* luaSt)
{
    static Spherical s = EMPTY_SPHERICAL;
    static const MWNamedArg argTable[] =
        {
            { "mass",  LUA_TNUMBER,  NULL, TRUE, &s.mass  },
            { "scale", LUA_TNUMBER,  NULL, TRUE, &s.scale },
            END_MW_NAMED_ARG
        };

    s.type = HernquistSpherical;
    return createSpherical(luaSt, argTable, &s);
}

static int createPlummerSpherical(lua_State* luaSt)
{
    static Spherical s = EMPTY_SPHERICAL;
    static const MWNamedArg argTable[] =
        {
            { "mass",  LUA_TNUMBER,  NULL, TRUE, &s.mass  },
            { "scale", LUA_TNUMBER,  NULL, TRUE, &s.scale },
            END_MW_NAMED_ARG
        };

    s.type = PlummerSpherical;
    return createSpherical(luaSt, argTable, &s);
}

static int createNoSpherical(lua_State* luaSt)
{
    static Spherical s = EMPTY_SPHERICAL;
    static const MWNamedArg argTable[] =
        {
            { "mass",  LUA_TNUMBER,  NULL, TRUE, &s.mass  },
            END_MW_NAMED_ARG
        };

    s.type = NoSpherical;
    return createSpherical(luaSt, argTable, &s);
}

int getSphericalT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, sphericalOptions, *(int*) v);
}

static int toStringSpherical(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showSpherical, (LuaTypeCheckFunc) checkSpherical);
}

static int eqSpherical(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalSpherical(checkSpherical(luaSt, 1), checkSpherical(luaSt, 2)));
    return 1;
}

int getSpherical(lua_State* luaSt, void* v)
{
    pushSpherical(luaSt, (Spherical*) v);
    return 1;
}

int setSpherical(lua_State* luaSt, void* v)
{
    *(Spherical*) v = *checkSpherical(luaSt, 3);
    return 0;
}

static const luaL_reg metaMethodsSpherical[] =
{
    { "__tostring", toStringSpherical },
    { "__eq",       eqSpherical       },
    { NULL, NULL }
};

static const luaL_reg methodsSpherical[] =
{
    { "hernquist", createHernquistSpherical },
    { "plummer", createPlummerSpherical },
    { "none", createNoSpherical },
    { NULL, NULL }
};

static const Xet_reg_pre gettersSpherical[] =
{
    { "type",  getSphericalT, offsetof(Spherical, type) },
    { "mass",  getNumber,     offsetof(Spherical, mass) },
    { "scale", getNumber,     offsetof(Spherical, scale) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersSpherical[] =
{
    { "mass",  setNumber, offsetof(Spherical, mass) },
    { "scale", setNumber, offsetof(Spherical, scale) },
    { NULL, NULL, 0 }
};

int registerSpherical(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          SPHERICAL_TYPE,
                          gettersSpherical,
                          settersSpherical,
                          metaMethodsSpherical,
                          methodsSpherical);
}

int registerSphericalKinds(lua_State* luaSt)
{
    int table;

    lua_newtable(luaSt);
    table = lua_gettop(luaSt);

    setModelTableItem(luaSt, table, createHernquistSpherical, "hernquist");
    setModelTableItem(luaSt, table, createPlummerSpherical, "plummer");
    setModelTableItem(luaSt, table, createNoSpherical, "none");

    lua_setglobal(luaSt, "sphericalModels");

    return 0;
}

