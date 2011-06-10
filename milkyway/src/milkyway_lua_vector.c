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

#include "milkyway_lua_marshal.h"
#include "milkyway_lua_vector.h"
#include "milkyway_util.h"


mwvector* checkVector(lua_State* luaSt, int idx)
{
    return (mwvector*) mw_checknamedudata(luaSt, idx, MWVECTOR_TYPE);
}

int pushVector(lua_State* luaSt, mwvector vIn)
{
    mwvector* vNew;

    vNew = (mwvector*) lua_newuserdata(luaSt, sizeof(mwvector));
    *vNew = vIn;

    luaL_getmetatable(luaSt, MWVECTOR_TYPE);
    lua_setmetatable(luaSt, -2);

    return 1;
}

int getVector(lua_State* luaSt, void* v)
{
    pushVector(luaSt, *(mwvector*) v);
    return 1;
}

int setVector(lua_State* luaSt, void* v)
{
    *(mwvector*)v = *checkVector(luaSt, 3);
    return 0;
}

static int absVector(lua_State* luaSt)
{
    mwvector* v;

    v = checkVector(luaSt, 1);
    lua_pushnumber(luaSt, (lua_Number) mw_absv(*v));

    return 1;
}

static int lengthVector(lua_State* luaSt)
{
    lua_pushnumber(luaSt, (lua_Number) mw_length(*checkVector(luaSt, 1)));
    return 1;
}

static const mwvector _zeroVector = ZERO_VECTOR;

static int createVector(lua_State* luaSt)
{
    int n;
    mwvector v;

    n = lua_gettop(luaSt);  /* Number of arguments */

    switch (n)
    {
        case 0:
            pushVector(luaSt, _zeroVector);
            return 1;

        case 3:
            v.x = luaL_checknumber(luaSt, 1);
            v.y = luaL_checknumber(luaSt, 2);
            v.z = luaL_checknumber(luaSt, 3);
            v.w = 0.0;
            pushVector(luaSt, v);
            return 1;

        default:
            return luaL_argerror(luaSt, 3, "Expected 0 or 3 arguments to create vector");
    }

    return 1;
}

static int toStringVector(lua_State* luaSt)
{
    mwvector* v;
    char* str;

    v = checkVector(luaSt, 1);
    str = showVector(*v);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

/* Get a scalar and a vector argument in either order for mixed
   vector/scalar maths.

   Returns  0 if not a vector and a scalar.
   Returns -1 if order is vector, scalar
   Returns +1 if order is scalar, vector
*/
static int checkScalarVectorArgs(lua_State* luaSt, real* s, mwvector* v)
{
    if (lua_isnumber(luaSt, 1))
    {
        *s = lua_tonumber(luaSt, 1);
        *v = *checkVector(luaSt, 2);
        return 1;
    }
    else if (lua_isnumber(luaSt, 2))
    {
        *v = *checkVector(luaSt, 1);
        *s = lua_tonumber(luaSt, 2);
        return -1;
    }
    else
    {
        return 0;
    }
}

static inline void check2Vector(lua_State* luaSt, mwvector* v1, mwvector* v2)
{
    *v1 = *checkVector(luaSt, 1);
    *v2 = *checkVector(luaSt, 2);
}

static int addVector(lua_State* luaSt)
{
    mwvector v1, v2;

    check2Vector(luaSt, &v1, &v2);
    pushVector(luaSt, mw_addv(v1, v2));
    return 1;
}

static int crossVector(lua_State* luaSt)
{
    mwvector v1, v2;

    check2Vector(luaSt, &v1, &v2);
    pushVector(luaSt, mw_crossv(v1, v2));
    return 1;
}

static int subVector(lua_State* luaSt)
{
    mwvector v1, v2;

    check2Vector(luaSt, &v1, &v2);
    pushVector(luaSt, mw_subv(v1, v2));
    return 1;
}

static int distVector(lua_State* luaSt)
{
    mwvector v1, v2;

    check2Vector(luaSt, &v1, &v2);
    lua_pushnumber(luaSt, mw_distv(v1, v2));
    return 1;
}

static int divVector(lua_State* luaSt)
{
    mwvector v1, v2;
    real s;

    /* What would it mean to divide a scalar by a vector? */
    if (lua_isnumber(luaSt, 2))
    {
        v1 = *checkVector(luaSt, 1);
        s = lua_tonumber(luaSt, 2);
        pushVector(luaSt, mw_divvs(v1, s));
    }
    else
    {
        check2Vector(luaSt, &v1, &v2);
        pushVector(luaSt, mw_divv(v1, v2));
    }

    return 1;
}

static int multVector(lua_State* luaSt)
{
    mwvector v1, v2;
    real s;
    int rc;

    rc = checkScalarVectorArgs(luaSt, &s, &v1);
    if (rc == 0)
    {
        /* Vector * Vector */
        check2Vector(luaSt, &v1, &v2);
        lua_pushnumber(luaSt, mw_dotv(v1, v2));
    }
    else
    {
        /* Vector * Scalar */
        pushVector(luaSt, mw_mulvs(v1, s));
    }

    return 1;
}

static int negVector(lua_State* luaSt)
{
    mwvector v;

    v = *checkVector(luaSt, 1);
    pushVector(luaSt, mw_negv(v));

    return 1;
}

static const luaL_reg metaMethodsVector[] =
{
    { "__tostring", toStringVector },
    { "__add",      addVector      },
    { "__sub",      subVector      },
    { "__mul",      multVector     },
    { "__div",      divVector      },
    { "__unm",      negVector      },
    { NULL, NULL }
};

static const luaL_reg methodsVector[] =
{
    { "create", createVector },
    { "abs",    absVector},
    { "length", lengthVector},
    { "cross",  crossVector},
    { "dist",   distVector},
    { NULL, NULL }
};

static const Xet_reg_pre gettersVector[] =
{
    { "x", getNumber, offsetof(mwvector, x) },
    { "y", getNumber, offsetof(mwvector, y) },
    { "z", getNumber, offsetof(mwvector, z) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersVector[] =
{
    { "x", setNumber, offsetof(mwvector, x) },
    { "y", setNumber, offsetof(mwvector, y) },
    { "z", setNumber, offsetof(mwvector, z) },
    { NULL, NULL, 0 }
};


int registerVector(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          MWVECTOR_TYPE,
                          gettersVector,
                          settersVector,
                          metaMethodsVector,
                          methodsVector);
}

