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
#include "nbody_lua_types.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"
#include "nbody_show.h"

static int totalBodies(lua_State* luaSt, int nModels)
{
    int top, i, n = 0;

    top = lua_gettop(luaSt);
    for (i = top; i > top - nModels; --i)
    {
        if (expectTable(luaSt, i))
        {
            mw_lua_perror(luaSt, "Error reading body table");
            return 0;
        }

        n += luaL_getn(luaSt, i);
    }

    return n;
}

static int readBodyArray(lua_State* luaSt, int table, Body* bodies, int n)
{
    int i;
    Body* b;

    for (i = 0; i < n; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);
        b = expectBody(luaSt, lua_gettop(luaSt));
        if (!b)
        {
            mw_lua_perror(luaSt, "Error reading body %d", i);
            break;
        }

        bodies[i] = *b;
        lua_pop(luaSt, 1);
    }

    return i != n; /* Didn't read all bodies successfully */
}

/* Read returned table of model components. Pops the n arguments */
Body* readModels(lua_State* luaSt, int nModels, int* nOut)
{
    int i, n, totalN, top;
    Body* allBodies;
    Body* bodies;

    totalN = totalBodies(luaSt, nModels);
    if (totalN == 0)
    {
        mw_printf("Didn't get any bodies\n");
        return NULL;
    }

    bodies = allBodies = (Body*) mwCallocA(totalN, sizeof(Body));

    for (i = 0; i < nModels; ++i)
    {
        top = lua_gettop(luaSt);
        n = luaL_getn(luaSt, top);

        if (readBodyArray(luaSt, top, bodies, n))
        {
            mw_printf("Error reading body array %d\n", i);
            free(allBodies);
            allBodies = NULL;
            totalN = 0;
            break;
        }

        bodies = &bodies[n];
        lua_pop(luaSt, 1);
    }

    lua_pop(luaSt, nModels - i);  /* No more body tables if error */

    if (nOut)
        *nOut = totalN;

    return allBodies;
}

