
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

#define NBODY_CTX "NBodyCtx"

NBodyCtx* checkNBodyCtx(lua_State* luaSt, int index)
{
    NBodyCtx* ctx;

    luaL_checktype(luaSt, index, LUA_TUSERDATA);

    ctx = (NBodyCtx*) luaL_checkudata(luaSt, index, NBODY_CTX);
    if (!ctx)
        luaL_typerror(luaSt, index, NBODY_CTX);

    return ctx;
}

static int positionNBodyCtx(lua_State *L)
{
    NBodyCtx* ctx;

    ctx = checkNBodyCtx(L, 1);

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

NBodyCtx* pushNBodyCtx(lua_State* luaSt)
{
    NBodyCtx* ctx;

    ctx = (NBodyCtx*)lua_newuserdata(luaSt, sizeof(NBodyCtx));
    luaL_getmetatable(luaSt, NBODY_CTX);
    lua_setmetatable(luaSt, -2);

    return ctx;
}

static int createNBodyCtx(lua_State* luaSt)
{
    NBodyCtx* ctx;

    ctx = pushNBodyCtx(luaSt);

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
    {"__gc", destroyNBodyCtx },
    {0, 0}
};

static const luaL_reg methodsNBodyCtx[] =
{
    {"create",   createNBodyCtx },
//    {"position", your_position},
    {"test",     test},
    {0, 0}
};

static const Xet_reg_pre gettersNBodyCtx[] =
{
    {"nbody",          getInt,    offsetof(NBodyCtx, nbody)          },
    {"timestep",       getNumber, offsetof(NBodyCtx, timestep)       },
    {"orbit_timestep", getNumber, offsetof(NBodyCtx, orbit_timestep) },
    {"sunGCDist",      getNumber, offsetof(NBodyCtx, sunGCDist)      },
    {0, 0, 0 }
};

static const Xet_reg_pre settersNBodyCtx[] =
{
    {"nbody",          setInt,    offsetof(NBodyCtx, nbody)          },
    {"timestep",       setNumber, offsetof(NBodyCtx, timestep)       },
    {"orbit_timestep", setNumber, offsetof(NBodyCtx, orbit_timestep) },
    {"sunGCDist",      setNumber, offsetof(NBodyCtx, sunGCDist)      },
    {0, 0, 0 }
};


int registerNBodyCtx(lua_State* luaSt)
{
    int metatable, methods;

    /* create methods table, & add it to the table of globals */
    luaL_openlib(luaSt, NBODY_CTX, methodsNBodyCtx, 0);
    methods = lua_gettop(luaSt);

    /* create metatable for NBodyCtx, & add it to the registry */
    luaL_newmetatable(luaSt, NBODY_CTX);
    luaL_openlib(luaSt, 0, metaMethodsNBodyCtx, 0);  /* fill metatable */
    metatable = lua_gettop(luaSt);

    lua_pushliteral(luaSt, "__metatable");
    lua_pushvalue(luaSt, methods);    /* dup methods table*/
    lua_rawset(luaSt, metatable);     /* hide metatable:
                                         metatable.__metatable = methods */
    lua_pushliteral(luaSt, "__index");
    lua_pushvalue(luaSt, metatable);     /* upvalue index 1 */
    Xet_add(luaSt, gettersNBodyCtx);     /* fill metatable with getters */
    lua_pushvalue(luaSt, methods);       /* upvalue index 2 */
    lua_pushcclosure(luaSt, indexHandler, 2);
    lua_rawset(luaSt, metatable);        /* metatable.__index = indexHandler */

    lua_pushliteral(luaSt, "__newindex");
    lua_newtable(luaSt);                 /* table for members you can set */
    Xet_add(luaSt, settersNBodyCtx);     /* fill with setters */
    lua_pushcclosure(luaSt, newIndexHandler, 1);
    lua_rawset(luaSt, metatable);     /* metatable.__newindex = newIndexHandler */

    lua_pop(luaSt, 1);            /* drop metatable */
    return 1;                     /* return methods on the stack */
}

