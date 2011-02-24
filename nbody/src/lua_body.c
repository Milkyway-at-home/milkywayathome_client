
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "io.h"
#include "lua_type_marshal.h"
#include "lua_nbodyctx.h"
#include "lua_body.h"
#include "lua_vector.h"

#include "milkyway_util.h"

body* checkBody(lua_State* luaSt, int index)
{
    body* b;

    b = (body*) luaL_checkudata(luaSt, index, BODY_TYPE);
    luaL_argcheck(luaSt, b != NULL, 1, "`Body' expected");

    return b;
}

static int positionBody(lua_State *L, int index)
{
    body* ctx;

    ctx = checkBody(L, index);

    #if 0
    double   x = yd->x;
    double   y = yd->y;
    if (lua_gettop(L) > 1)
    {
        yd->x = luaL_checknumber(L, 2);
        yd->y = luaL_checknumber(L, 3);
    }
    lua_pushnumber(L,x);
    lua_pushnumber(L,y);

    #endif
    return 2;
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

static int destroyBody(lua_State* luaSt)
{
    body* b;

    b = (body*) lua_touserdata(luaSt, 1);
    printf("Goodbye body\n");
    return 0;
}

static const luaL_reg metaMethodsBody[] =
{
    { "__gc", destroyBody },
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

