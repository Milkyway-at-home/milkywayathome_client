
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

mwvector checkVector(lua_State* luaSt, int index)
{
    mwvector* v;

    v = (mwvector*) luaL_checkudata(luaSt, index, BODY_TYPE);
    luaL_argcheck(luaSt, v != NULL, 1, "`Vector' expected");

    return *v;
}

int pushVector(lua_State* luaSt, mwvector vIn)
{
    mwvector* vNew;

    vNew = (mwvector*) lua_newuserdata(luaSt, sizeof(mwvector));
    *vNew = vIn;

    //luaL_getmetatable(luaSt, MWVECTOR);
    //lua_setmetatable(luaSt, -2);

    return 1;
}

