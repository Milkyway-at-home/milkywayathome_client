
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

#include "milkyway_util.h"

body* checkBody(lua_State* luaSt, int index)
{
    body* b;

    b = (body*) luaL_checkudata(luaSt, index, BODY_TYPE);
    luaL_argcheck(luaSt, b != NULL, 1, "`Body' expected");

    return b;
}

int pushBody(lua_State* luaSt, const body* b)
{
    body* bArr;

    bArr = (body*) lua_newuserdata(luaSt, sizeof(body));
    *bArr = *b;

    //luaL_getmetatable(luaSt, BODY);
    //lua_setmetatable(luaSt, -2);

    return 1;
}

