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
#include "lua_body.h"
#include "lua_vector.h"

#include "milkyway_util.h"

body* checkBody(lua_State* luaSt, int index)
{
    return (body*) mw_checknamedudata(luaSt, index, BODY_TYPE);
}

int pushBody(lua_State* luaSt, const body* ctx)
{
    body* lctx;

    lctx = (body*)lua_newuserdata(luaSt, sizeof(body));
    if (!lctx)
    {
        warn("Creating Body userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, BODY_TYPE);
    lua_setmetatable(luaSt, -2);

    *lctx = *ctx;

    return 0;
}

static const body _emptyBody = EMPTY_BODY;


static int createBody(lua_State* luaSt)
{
    pushBody(luaSt, &_emptyBody);

    /* TODO: (mass, position, velocity) constructor */

   #if 0
    size_t name_len;
    const char* name;

    name = luaL_checklstring(luaSt, 1, &name_len);
    if (name_len > 15)
        luaL_error(luaSt, "name too long");

    strcpy(yd->name, name);
    yd->age = luaL_checkint(luaSt, 2);
    yd->x   = luaL_checknumber(luaSt, 3);
    yd->y   = luaL_checknumber(luaSt, 4);
    yd->id  = ++id_counter;
    #endif

    return 1;
}

static int toStringBody(lua_State* luaSt)
{
    body* b;
    char* str;

    b = checkBody(luaSt, 1);
    str = showBody(b);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static const luaL_reg metaMethodsBody[] =
{
    { "__tostring", toStringBody },
    { NULL, NULL }
};

static const luaL_reg methodsBody[] =
{
    { "create",   createBody },
    { NULL, NULL }
};

static const Xet_reg_pre gettersBody[] =
{
    { "velocity", getVector, offsetof(body, vel) },
    { "position", getVector, offsetof(body, bodynode.pos) },
    { "mass",     getNumber, offsetof(body, bodynode.mass) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersBody[] =
{
    { "velocity", setVector, offsetof(body, vel) },
    { "position", setVector, offsetof(body, bodynode.pos) },
    { "mass",     setNumber, offsetof(body, bodynode.mass) },
    { NULL, NULL, 0 }
};


int registerBody(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          BODY_TYPE,
                          gettersBody,
                          settersBody,
                          metaMethodsBody,
                          methodsBody);
}

