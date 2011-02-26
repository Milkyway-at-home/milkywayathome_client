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
#include "lua_vector.h"
#include "lua_initial_conditions.h"
#include "lua_nbodyctx.h"

#include "milkyway_util.h"
#include "nbody_util.h"


InitialConditions* checkInitialConditions(lua_State* luaSt, int index)
{
    return (InitialConditions*) mw_checknamedudata(luaSt, index, INITIAL_CONDITIONS_TYPE);
}

int pushInitialConditions(lua_State* luaSt, const InitialConditions* ic)
{
    InitialConditions* lic;

    lic = (InitialConditions*)lua_newuserdata(luaSt, sizeof(InitialConditions));
    if (!lic)
    {
        warn("Creating InitialConditions userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, INITIAL_CONDITIONS_TYPE);
    lua_setmetatable(luaSt, -2);

    *lic = *ic;

    return 0;
}

static int createInitialConditions(lua_State* luaSt)
{
    InitialConditions ic = EMPTY_INITIAL_CONDITIONS;

    warn("Creating initial conditions\n");
    pushInitialConditions(luaSt, &ic);
    return 1;
}

static int toStringInitialConditions(lua_State* luaSt)
{
    InitialConditions* ic;
    char* str;

    ic = checkInitialConditions(luaSt, 1);
    str = showInitialConditions(ic);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static const luaL_reg metaMethodsInitialConditions[] =
{
    { "__tostring", toStringInitialConditions },
    { NULL, NULL }
};

static const luaL_reg methodsInitialConditions[] =
{
    { "create",   createInitialConditions },
    { NULL, NULL }
};

/* TODO Error when writing to fields a halo type doesn't have */
static const Xet_reg_pre gettersInitialConditions[] =
{
    { "position",     getVector, offsetof(InitialConditions, position)     },
    { "velocity",     getVector, offsetof(InitialConditions, velocity)     },
    { "useGalC",      getBool,   offsetof(InitialConditions, useGalC)      },
    { "useRadians",   getBool,   offsetof(InitialConditions, useRadians)   },
    { "reverseOrbit", getBool,   offsetof(InitialConditions, reverseOrbit) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersInitialConditions[] =
{
    { "position",     setVector, offsetof(InitialConditions, position)     },
    { "velocity",     setVector, offsetof(InitialConditions, velocity)     },
    { "useGalC",      setBool,   offsetof(InitialConditions, useGalC)      },
    { "useRadians",   setBool,   offsetof(InitialConditions, useRadians)   },
    { "reverseOrbit", setBool,   offsetof(InitialConditions, reverseOrbit) },
    { NULL, NULL, 0 }
};

int registerInitialConditions(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          INITIAL_CONDITIONS_TYPE,
                          gettersInitialConditions,
                          settersInitialConditions,
                          metaMethodsInitialConditions,
                          methodsInitialConditions);
}

