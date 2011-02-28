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
#include "show.h"
#include "lua_type_marshal.h"
#include "lua_spherical.h"

#include "milkyway_util.h"

Spherical* checkSpherical(lua_State* luaSt, int index)
{
    return (Spherical*) mw_checknamedudata(luaSt, index, SPHERICAL_TYPE);
}

int pushSpherical(lua_State* luaSt, const Spherical* d)
{
    Spherical* ld;

    ld = (Spherical*) lua_newuserdata(luaSt, sizeof(Spherical));
    if (!ld)
    {
        warn("Creating Spherical userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, SPHERICAL_TYPE);
    lua_setmetatable(luaSt, -2);

    *ld = *d;

    return 0;
}

static const MWEnumAssociation sphericalOptions[] =
{
    { "spherical", SphericalPotential },
    END_MW_ENUM_ASSOCIATION
};

static int createSpherical(lua_State* luaSt)
{
    Spherical s;
    int nArgs;
    real mass = NAN, scale = NAN;
    spherical_t type = SphericalPotential;

    const MWNamedArg argTable[] =
        {
          //{ "type",  ,  NULL, FALSE, &type  },
            { "mass",  LUA_TNUMBER,  NULL, FALSE, &mass  },
            { "scale", LUA_TNUMBER,  NULL, FALSE, &scale },
            END_MW_NAMED_ARG
        };

    /* Enum option isn't useful now */
    nArgs = lua_gettop(luaSt);
    if (nArgs == 1 && lua_istable(luaSt, 1))
    {
        handleNamedArgumentTable(luaSt, argTable, 1);
    }
    else if (nArgs == 1 || nArgs == 2)
    {
        mass = luaL_optnumber(luaSt, 1, mass);
        scale = luaL_optnumber(luaSt, 2, scale);
    }
    else
    {
        return luaL_argerror(luaSt, 1, "Expected 1 or 2 arguments");
    }

    s.type = type;
    s.mass = mass;
    s.scale = scale;
    pushSpherical(luaSt, &s);

    return 1;
}

int getSphericalT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, sphericalOptions, *(int*) v);
}

static int toStringSpherical(lua_State* luaSt)
{
    Spherical* d;
    char* str;

    d = checkSpherical(luaSt, 1);
    str = showSpherical(d);
    lua_pushstring(luaSt, str);
    free(str);

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
    { NULL, NULL }
};

static const luaL_reg methodsSpherical[] =
{
    { "create", createSpherical },
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
    //{ "type",  setSphericalT, offsetof(Spherical, type) },
    { "mass",  setNumber,     offsetof(Spherical, mass) },
    { "scale", setNumber,     offsetof(Spherical, scale) },
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

