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
#include "nbody_potential.h"
#include "nbody_lua_potential.h"
#include "nbody_lua_disk.h"
#include "nbody_lua_halo.h"
#include "nbody_lua_spherical.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"

Potential* checkPotential(lua_State* luaSt, int idx)
{
    return (Potential*) mw_checknamedudata(luaSt, idx, POTENTIAL_TYPE);
}

int pushPotential(lua_State* luaSt, const Potential* p)
{
    return pushType(luaSt, POTENTIAL_TYPE, sizeof(Potential), (void*) p);
}

Potential* toPotential(lua_State* luaSt, int idx)
{
    return (Potential*) mw_tonamedudata(luaSt, idx, POTENTIAL_TYPE);
}

Potential* expectPotential(lua_State* luaSt, int idx)
{
    return (Potential*) expectType(luaSt, idx, POTENTIAL_TYPE);
}

int getPotential(lua_State* luaSt, void* v)
{
    return pushPotential(luaSt, (Potential*) v);
}

int setPotential(lua_State* luaSt, void* v)
{
    *(Potential*) v = *checkPotential(luaSt, 3);
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

static int luaAcceleration(lua_State* luaSt)
{
    const Potential* pot;
    const mwvector* r;

    pot = checkPotential(luaSt, 1);
    r = checkVector(luaSt, 2);

    pushVector(luaSt, nbExtAcceleration(pot, *r));
    return 1;
}

static int toStringPotential(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showPotential, (LuaTypeCheckFunc) checkPotential);
}

static const luaL_reg metaMethodsPotential[] =
{
    { "__tostring", toStringPotential },
    { NULL, NULL }
};

static const luaL_reg methodsPotential[] =
{
    { "create",       createPotential },
    { "acceleration", luaAcceleration },
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

