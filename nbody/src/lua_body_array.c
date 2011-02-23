
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
#include "lua_body_array.h"

#include "milkyway_util.h"


NBodyLuaBodyArray* checkNBodyLuaBodyArray(lua_State* luaSt, int index)
{
    void* ud;

    ud = (NBodyLuaBodyArray*) luaL_checkudata(luaSt, index, NBODY_LUA_BODY_ARRAY);
    luaL_argcheck(luaSt, ud != NULL, 1, "`NBodyLuaBodyArray' expected");

    return ud;
}

static int getSizeNBodyLuaBodyArray(lua_State* luaSt)
{
    NBodyLuaBodyArray* arr;

    warn("Size\n");

    arr = checkNBodyLuaBodyArray(luaSt, -1);
    lua_pushnumber(luaSt, arr->nBody);

    return 1;
}

static body* getBodyNBodyLuaBodyArray(lua_State* luaSt)
{
    int index;
    NBodyLuaBodyArray* arr;

    warn("Element Getter\n");

    arr = checkNBodyLuaBodyArray(luaSt, 1);
    index = luaL_checkint(luaSt, 2);
    luaL_argcheck(luaSt, 1 <= index && index <= arr->nBody, 2, "index out of range");

    warn("Got some elements\n");

    return &arr->bodies[index - 1];
}

static int indexHandlerNBodyLuaBodyArray(lua_State* luaSt)
{
    /* stack has userdata, index */
    int idx;
    NBodyLuaBodyArray* arr;

    warn("soup\n");
    idx = luaL_checkint(luaSt, -1);
    //arr = checkNBodyLuaBodyArray(luaSt, -2);

    if (idx > arr->nBody || idx < 0)
        luaL_error(luaSt, "Body index (%d) out of bounds (%d)\n", idx, arr->nBody);

    // push body

    return 0;
}

static int getNBodyLuaBodyArray(lua_State* luaSt)
{
    body* b;

    b = getBodyNBodyLuaBodyArray(luaSt);
    pushBody(luaSt, b);

    return 1;
}

static int setNBodyLuaBodyArray(lua_State* luaSt)
{
    body* newB;
    body* oldB;

    newB = checkBody(luaSt, 3);
    oldB = getBodyNBodyLuaBodyArray(luaSt);

    *oldB = *newB;

    return 0;
}


static int newNBodyLuaBodyArray(lua_State* luaSt)
{
    int n;
    NBodyLuaBodyArray* arr;

    n = luaL_checkint(luaSt, 1);
    arr = (NBodyLuaBodyArray*) lua_newuserdata(luaSt, sizeof(NBodyLuaBodyArray));

    luaL_getmetatable(luaSt, NBODY_LUA_BODY_ARRAY);
    lua_setmetatable(luaSt, -2);

    arr->nBody = 0;
    arr->bodies = NULL;

    if (n <= 0)
        luaL_error(luaSt, "Number of bodies must not be <= 0");
    arr->nBody = n;
    arr->bodies = mwCalloc(arr->nBody, sizeof(body));

    return 1;
}

static int stringFromNBodyLuaBodyArray(lua_State* luaSt)
{
    NBodyLuaBodyArray* arr;

    arr = checkNBodyLuaBodyArray(luaSt, -1);
    lua_pushfstring(luaSt, "NBodyLuaBodyArray(%d)", arr->nBody);
    return 1;
}

static const struct luaL_reg methodsNBodyLuaBodyArray [] =
{
    { "new", newNBodyLuaBodyArray },
    { NULL, NULL }
};

static const struct luaL_reg instanceMethodsNBodyLuaBodyArray [] =
{
    { "__tostring", stringFromNBodyLuaBodyArray },
    { "set",        setNBodyLuaBodyArray        },
    { "get",        getNBodyLuaBodyArray        },
    { "size",       getSizeNBodyLuaBodyArray    },
    { NULL, NULL }
};

static int destroyNBodyLuaBodyArray(lua_State* luaSt)
{
    return 0;
}

static const luaL_reg metaMethodsNBodyLuaBodyArray[] =
{
    { "__gc", destroyNBodyLuaBodyArray },
    { 0, 0 }
};


static void setNBodyLuaBodyArrayIndexHandlers(lua_State* luaSt, int metatable, int methods)
{
    lua_pushliteral(luaSt, "__index");
    lua_pushvalue(luaSt, metatable);     /* upvalue index 1 */
    lua_pushvalue(luaSt, methods);       /* upvalue index 2 */
    lua_pushcclosure(luaSt, indexHandlerNBodyLuaBodyArray, 2);
    lua_rawset(luaSt, metatable);        /* metatable.__index = indexHandler */

    #if 0
    lua_pushliteral(luaSt, "__newindex");
    //lua_newtable(luaSt);                 /* table for members you can set */
    //Xet_add(luaSt, settersNBodyCtx);     /* fill with setters */
    lua_pushcclosure(luaSt, newIndexHandler, 1);
    lua_rawset(luaSt, metatable);     /* metatable.__newindex = newIndexHandler */
    #endif
}

static int openNBodyLuaBodyArray(lua_State* luaSt)
{
    luaL_newmetatable(luaSt, NBODY_LUA_BODY_ARRAY);
    lua_pushstring(luaSt, "__index");
    lua_pushvalue(luaSt, -2);   /* Pushes the metatable */
    lua_settable(luaSt, -3);    /* metatable.__index = metatable */

    luaL_openlib(luaSt, NULL, instanceMethodsNBodyLuaBodyArray, 0);
    luaL_openlib(luaSt, NBODY_LUA_BODY_ARRAY_LIB, methodsNBodyLuaBodyArray, 0);

    return 1;
}

int registerNBodyLuaBodyArray(lua_State* luaSt)
{
    //int metatable, methods;

    warn("Register NBodyLuaBodyArray\n");

    return openNBodyLuaBodyArray(luaSt);

    //lua_pushliteral(luaSt, "__metatable");
    //lua_pushvalue(luaSt, methods);    /* dup methods table*/
    #if 0
    lua_rawset(luaSt, metatable);     /* hide metatable:
                                         metatable.__metatable = methods */
#endif
    //setNBodyLuaBodyArrayIndexHandlers(luaSt, metatable, methods);

    //lua_pop(luaSt, 1);            /* drop metatable */
    return 1;                     /* return methods on the stack */
}



