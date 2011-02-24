
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "io.h"
#include "lua_type_marshal.h"
#include "lua_vector.h"

#include "milkyway_util.h"

mwvector* checkVector(lua_State* luaSt, int index)
{
    mwvector* v;

    v = (mwvector*) luaL_checkudata(luaSt, index, MWVECTOR);
    luaL_argcheck(luaSt, v != NULL, 1, "`Vector' expected");

    return v;
}

int pushVector(lua_State* luaSt, mwvector vIn)
{
    mwvector* vNew;

    vNew = (mwvector*) lua_newuserdata(luaSt, sizeof(mwvector));
    *vNew = vIn;

    luaL_getmetatable(luaSt, MWVECTOR);
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

    v = checkVector(luaSt, -1);
    lua_pushnumber(luaSt, (lua_Number) mw_absv(*v));

    return 1;
}

static const mwvector _emptyVector = EMPTY_MWVECTOR;

static int createVector(lua_State* luaSt)
{
    int n;
    mwvector v;

    n = lua_gettop(luaSt);  /* Number of arguments */

    switch (n)
    {
        case 0:
            warn("Creating 0 arg vector\n");
            pushVector(luaSt, _emptyVector);
            return 1;

        case 3:
            warn("Creating 3 arg vector\n");
            v.x = luaL_checknumber(luaSt, 1);
            v.y = luaL_checknumber(luaSt, 2);
            v.z = luaL_checknumber(luaSt, 3);
            v.w = 0.0;
            pushVector(luaSt, v);
            return 1;

        default:
            warn("Vector constructor failure\n");
            luaL_argerror(luaSt, 3, "Expected 0 or 3 arguments to create vector");
            return 0;
    }

    if (n == 0)
    {
    }

    // pushVector(luaSt, _emptyVector);

    return 3;
}

static int destroyVector(lua_State* luaSt)
{
    mwvector* ctx;

    ctx = (mwvector*) lua_touserdata(luaSt, 1);
    printf("Goodbye vector\n");
    return 0;
}

static const luaL_reg metaMethodsVector[] =
{
    { NULL, NULL }
};

static const luaL_reg methodsVector[] =
{
    { "create", createVector },
    { "abs",    absVector},
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
                          MWVECTOR,
                          gettersVector,
                          settersVector,
                          metaMethodsVector,
                          methodsVector);
}

