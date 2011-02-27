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
#include "lua_halo.h"

#include "milkyway_util.h"

Halo* checkHalo(lua_State* luaSt, int index)
{
    return (Halo*) mw_checknamedudata(luaSt, index, HALO_TYPE);
}

int pushHalo(lua_State* luaSt, const Halo* h)
{
    Halo* lh;

    lh = (Halo*)lua_newuserdata(luaSt, sizeof(Halo));
    if (!lh)
    {
        warn("Creating Halo userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, HALO_TYPE);
    lua_setmetatable(luaSt, -2);

    *lh = *h;

    return 0;
}

static const MWEnumAssociation haloOptions[] =
{
    { "logarithmic", LogarithmicHalo },
    { "nfw",         NFWHalo,        },
    { "triaxial",    TriaxialHalo,   },
    END_MW_ENUM_ASSOCIATION
};

static int createHalo(lua_State* luaSt)
{
    Halo h = EMPTY_HALO;

    warn("Creating halo\n");

    h.type = checkEnum(luaSt, haloOptions, -1);
    pushHalo(luaSt, &h);
    return 1;
}

static int toStringHalo(lua_State* luaSt)
{
    Halo* h;
    char* str;

    h = checkHalo(luaSt, 1);
    str = showHalo(h);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

int getHaloT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, haloOptions, *(int*) v);
}

static const luaL_reg metaMethodsHalo[] =
{
    { "__tostring", toStringHalo },
    { NULL, NULL }
};

static const luaL_reg methodsHalo[] =
{
    { "create",   createHalo },
    { NULL, NULL }
};

/* TODO Error when writing to fields a halo type doesn't have */
static const Xet_reg_pre gettersHalo[] =
{
    { "type",         getHaloT,  offsetof(Halo, type) },
    { "vhalo",        getNumber, offsetof(Halo, vhalo) },
    { "scale_length", getNumber, offsetof(Halo, scale_length) },
    { "flattenX",     getNumber, offsetof(Halo, flattenX) },
    { "flattenY",     getNumber, offsetof(Halo, flattenY) },
    { "flattenZ",     getNumber, offsetof(Halo, flattenZ) },
    { "triaxAngle",   getNumber, offsetof(Halo, triaxAngle) },
    { "c1",           getNumber, offsetof(Halo, c1) },
    { "c2",           getNumber, offsetof(Halo, c2) },
    { "c3",           getNumber, offsetof(Halo, c3) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersHalo[] =
{
    { "vhalo",        setNumber, offsetof(Halo, vhalo) },
    { "scale_length", setNumber, offsetof(Halo, scale_length) },
    { "flattenX",     setNumber, offsetof(Halo, flattenX) },
    { "flattenY",     setNumber, offsetof(Halo, flattenY) },
    { "flattenZ",     setNumber, offsetof(Halo, flattenZ) },
    { "triaxAngle",   setNumber, offsetof(Halo, triaxAngle) },
    { "c1",           setNumber, offsetof(Halo, c1) },
    { "c2",           setNumber, offsetof(Halo, c2) },
    { "c3",           setNumber, offsetof(Halo, c3) },
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

