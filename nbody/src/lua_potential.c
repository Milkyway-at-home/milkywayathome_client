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
#include "lua_potential.h"
#include "lua_disk.h"
#include "lua_halo.h"
#include "lua_spherical.h"

#include "milkyway_util.h"

Potential* checkPotential(lua_State* luaSt, int index)
{
    return (Potential*) mw_checknamedudata(luaSt, index, POTENTIAL_TYPE);
}

int pushPotential(lua_State* luaSt, const Potential* p)
{
    Potential* lp;

    lp = (Potential*) lua_newuserdata(luaSt, sizeof(Potential));
    if (!lp)
    {
        warn("Creating Potential userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, POTENTIAL_TYPE);
    lua_setmetatable(luaSt, -2);

    *lp = *p;

    return 0;
}

static int createPotential(lua_State* luaSt)
{
    Potential p = EMPTY_POTENTIAL;
    const Spherical* s = NULL;
    const Disk* d = NULL;
    const Halo* h = NULL;

    const MWNamedArg argTable[] =
        {
            { "spherical", LUA_TUSERDATA, SPHERICAL_TYPE, TRUE,  &s },
            { "halo",      LUA_TUSERDATA, HALO_TYPE,      TRUE,  &h },
            { "disk",      LUA_TUSERDATA, DISK_TYPE,      TRUE,  &d },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:  /* Table of named arguments */
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 3:
            s = checkSpherical(luaSt, 1);
            d = checkDisk(luaSt, 2);
            h = checkHalo(luaSt, 3);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 3 arguments");
    }

    p.sphere[0] = *s;
    p.disk = *d;
    p.halo = *h;

    pushPotential(luaSt, &p);
    return 1;
}

static int toStringPotential(lua_State* luaSt)
{
    Potential* p;
    char* str;

    p = checkPotential(luaSt, 1);
    str = showPotential(p);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static const luaL_reg metaMethodsPotential[] =
{
    { "__tostring", toStringPotential },
    { NULL, NULL }
};

static const luaL_reg methodsPotential[] =
{
    { "create",   createPotential },
    { NULL, NULL }
};

static const Xet_reg_pre gettersPotential[] =
{
    { "spherical", getSpherical, offsetof(Potential, sphere[0]) },
    { "disk" ,     getDisk,      offsetof(Potential, disk)      },
    { "halo",      getHalo,      offsetof(Potential, halo)      },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersPotential[] =
{
    { "spherical", setSpherical, offsetof(Potential, sphere[0]) },
    { "disk" ,     setDisk,      offsetof(Potential, disk)      },
    { "halo",      setHalo,      offsetof(Potential, halo)      },
    { NULL, NULL, 0 }
};

int registerPotential(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          POTENTIAL_TYPE,
                          gettersPotential,
                          settersPotential,
                          metaMethodsPotential,
                          methodsPotential);
}

