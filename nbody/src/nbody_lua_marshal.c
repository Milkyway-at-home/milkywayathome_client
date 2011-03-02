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

#include "lua_type_marshal.h"
#include "nbody_lua_types.h"
#include "nbody_types.h"
#include "milkyway_util.h"
#include "nbody_lua_marshal.h"


static int totalBodies(lua_State* luaSt, int nModels)
{
    int top, i, n = 0;

    top = lua_gettop(luaSt);
    for (i = top; i > top - nModels; --i)
    {
        luaL_checktype(luaSt, i, LUA_TTABLE);
        n += luaL_getn(luaSt, i);
    }

    return n;
}

static void readBodyArray(lua_State* luaSt, int table, body* bodies, int n)
{
    int i;

    for (i = 0; i < n; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);
        bodies[i] = *checkBody(luaSt, -1);
        lua_pop(luaSt, 1);
    }
}

/* Read returned table of model components. Pops the n arguments */
body* readReturnedModels(lua_State* luaSt, int nModels, unsigned int* nOut)
{
    int i, n, totalN, top;
    body* allBodies;
    body* bodies;

    totalN = totalBodies(luaSt, nModels);
    if (totalN == 0)
    {
        warn("Didn't get any bodies\n");
        return NULL;
    }

    bodies = allBodies = (body*) mwMallocA(totalN * sizeof(body));

    for (i = 0; i < nModels; ++i)
    {
        top = lua_gettop(luaSt);
        n = luaL_getn(luaSt, top);

        readBodyArray(luaSt, top, bodies, n);
        bodies = &bodies[n];

        lua_pop(luaSt, 1);  /* No more body table */
    }

    if (nOut)
        *nOut = (unsigned int) totalN;

    return allBodies;
}


