
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

#include "milkyway_util.h"

#define NBODY_CTX "NBodyCtx"

NBodyCtx* checkNBodyCtx(lua_State* luaSt, int index)
{
    NBodyCtx* ctx;

    ctx = (NBodyCtx*) luaL_checkudata(luaSt, index, NBODY_CTX);
    luaL_argcheck(luaSt, ctx != NULL, 1, "`NBodyCtx' expected");

    return ctx;
}

static int positionNBodyCtx(lua_State *L, int index)
{
    NBodyCtx* ctx;

    ctx = checkNBodyCtx(L, index);

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

    return 2;

    #endif
    return 0;
}

int pushNBodyCtx(lua_State* luaSt, const NBodyCtx* ctx)
{
    NBodyCtx* lctx;

    lctx = (NBodyCtx*)lua_newuserdata(luaSt, sizeof(NBodyCtx));
    if (!lctx)
    {
        warn("Creating NBodyCtx userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, NBODY_CTX);
    lua_setmetatable(luaSt, -2);

    *lctx = *ctx;

    return 0;
}

static const NBodyCtx _emptyCtx = EMPTY_CTX;

static int createNBodyCtx(lua_State* luaSt)
{
    pushNBodyCtx(luaSt, &_emptyCtx);
    return 1;
}

static int destroyNBodyCtx(lua_State* luaSt)
{
    NBodyCtx* ctx;

    ctx = (NBodyCtx*) lua_touserdata(luaSt, 1);
    printf("Goodbye: %f\n", ctx->timestep);
    return 0;
}

static int test(lua_State* luaSt)
{
    int n;

    n = luaL_checknumber(luaSt, 1);
    lua_pushnumber(luaSt, 66);
    lua_pushnumber(luaSt, 67);
    lua_pushnumber(luaSt, 68);

    return n;
}

static const luaL_reg metaMethodsNBodyCtx[] =
{
    { "__gc", destroyNBodyCtx },
    { NULL, NULL }
};

static const luaL_reg methodsNBodyCtx[] =
{
    { "create",   createNBodyCtx },
//    {"position", your_position},
    { "test",     test },
    { NULL, NULL }
};

static const Xet_reg_pre gettersNBodyCtx[] =
{
    { "nbody",          getInt,    offsetof(NBodyCtx, nbody)          },
    { "timestep",       getNumber, offsetof(NBodyCtx, timestep)       },
    { "orbit_timestep", getNumber, offsetof(NBodyCtx, orbit_timestep) },
    { "sunGCDist",      getNumber, offsetof(NBodyCtx, sunGCDist)      },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersNBodyCtx[] =
{
    { "nbody",          setInt,    offsetof(NBodyCtx, nbody)          },
    { "timestep",       setNumber, offsetof(NBodyCtx, timestep)       },
    { "orbit_timestep", setNumber, offsetof(NBodyCtx, orbit_timestep) },
    { "sunGCDist",      setNumber, offsetof(NBodyCtx, sunGCDist)      },
    { NULL, NULL, 0 }
};


int registerNBodyCtx(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          NBODY_CTX,
                          gettersNBodyCtx,
                          settersNBodyCtx,
                          metaMethodsNBodyCtx,
                          methodsNBodyCtx);
}

