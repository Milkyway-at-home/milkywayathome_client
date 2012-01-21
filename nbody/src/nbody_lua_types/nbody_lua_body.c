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
#include "nbody_lua_body.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"

Body* toBody(lua_State* luaSt, int idx)
{
    return (Body*) mw_tonamedudata(luaSt, idx, BODY_TYPE);
}

Body* checkBody(lua_State* luaSt, int idx)
{
    return (Body*) mw_checknamedudata(luaSt, idx, BODY_TYPE);
}

/* Like checkBody() except for non-Lua called C functions */
Body* expectBody(lua_State* luaSt, int idx)
{
    Body* b;

    b = toBody(luaSt, idx);
    if (!b)
    {
        /* Push error message onto stack */
        mw_lua_typecheck(luaSt, idx, LUA_TUSERDATA, BODY_TYPE);
    }

    return b;
}


int pushBody(lua_State* luaSt, const Body* p)
{
    return pushType(luaSt, BODY_TYPE, sizeof(Body), (void*) p);
}

static const Body _emptyBody = EMPTY_BODY;


static int createBody(lua_State* luaSt)
{
    static Body b = EMPTY_BODY;
    static mwvector* x = NULL;
    static mwvector* v = NULL;
    static mwbool ignore = FALSE;
    static const MWNamedArg argTable[] =
        {
            { "mass",     LUA_TNUMBER,   NULL,          TRUE,  &b.bodynode.mass },
            { "position", LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &x               },
            { "velocity", LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &v               },
            { "ignore",   LUA_TBOOLEAN,  NULL,          FALSE, &ignore          },
            END_MW_NAMED_ARG
        };

    memset(&b, 0, sizeof(b));
    ignore = FALSE; /* Set to default value for table */
    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 3:
        case 4:
            b.bodynode.mass = luaL_checknumber(luaSt, 1);
            x = checkVector(luaSt, 2);
            v = checkVector(luaSt, 3);
            ignore = mw_lua_optboolean(luaSt, 4, FALSE);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1, 3 or 4 arguments");
    }

    b.bodynode.pos = *x;
    b.vel = *v;
    b.bodynode.type = BODY(ignore);

    pushBody(luaSt, &b);
    return 1;
}

static int getBodyIgnore(lua_State* luaSt, void* v)
{
    lua_pushboolean(luaSt, bodyTypeIsIgnore(*(body_t*) v));
    return 1;
}

static int setBodyIgnore(lua_State* luaSt, void* v)
{
    *(body_t*) v = BODY(mw_lua_checkboolean(luaSt, 3));
    return 0;
}

static int toStringBody(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showBody, (LuaTypeCheckFunc) checkBody);
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
    { "velocity", getVector,     offsetof(Body, vel)           },
    { "position", getVector,     offsetof(Body, bodynode.pos)  },
    { "mass",     getNumber,     offsetof(Body, bodynode.mass) },
    { "ignore",   getBodyIgnore, offsetof(Body, bodynode.type) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersBody[] =
{
    { "velocity", setVector,     offsetof(Body, vel)           },
    { "position", setVector,     offsetof(Body, bodynode.pos)  },
    { "mass",     setNumber,     offsetof(Body, bodynode.mass) },
    { "ignore",   setBodyIgnore, offsetof(Body, bodynode.type) },
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

