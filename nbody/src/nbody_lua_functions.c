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


#include "nbody_types.h"
#include "lua_type_marshal.h"
#include "milkyway_util.h"

#include "nbody_lua_functions.h"


static int luaGeneratePlummer(lua_State* luaSt)
{
    int nArgs;

    nArgs = lua_gettop(luaSt);

    warn("PLUMMER TEST nArgs = %d\n", nArgs);
    switch (nArgs)
    {
        case 1: /* Hopefully got table for fake named arguments */

        case 6: /* ignore argument is optional */

        case 7:

        default:
            warn("arst\n");
    }

    return 0;
}

static void registerLuaGeneratePlummer(lua_State* luaSt)
{
    lua_pushcfunction(luaSt, luaGeneratePlummer);
    lua_setglobal(luaSt, "generatePlummer");
}

void registerPredefinedModelGenerators(lua_State* luaSt)
{
    registerLuaGeneratePlummer(luaSt);
}


