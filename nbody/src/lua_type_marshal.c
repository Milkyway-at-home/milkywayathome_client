
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include "lua_type_marshal.h"
#include "nbody_types.h"

/* Mostly from example at http://lua-users.org/wiki/BindingWithMembersAndMethods */


int getInt(lua_State* luaSt, void* v)
{
    lua_pushnumber(luaSt, *(int*)v);
    return 1;
}

int setInt(lua_State* luaSt, void* v)
{
    *(int*)v = luaL_checkint(luaSt, 3);
    return 0;
}

int getNumber(lua_State* luaSt, void* v)
{
    lua_pushnumber(luaSt, *(lua_Number*)v);
    return 1;
}

int setNumber(lua_State* luaSt, void* v)
{
    *(lua_Number*)v = luaL_checknumber(luaSt, 3);
    return 0;
}

int getString(lua_State* luaSt, void* v)
{
    lua_pushstring(luaSt, (char*)v );
    return 1;
}

void Xet_add(lua_State* luaSt, Xet_reg l)
{
    for (; l->name; l++)
    {
        lua_pushstring(luaSt, l->name);
        lua_pushlightuserdata(luaSt, (void*)l);
        lua_settable(luaSt, -3);
    }
}

int Xet_call(lua_State* luaSt)
{
    /* for get: stack has userdata, index, lightuserdata */
    /* for set: stack has userdata, index, value, lightuserdata */
    Xet_reg m = (Xet_reg)lua_touserdata(luaSt, -1);  /* member info */
    lua_pop(luaSt, 1);                               /* drop lightuserdata */
    luaL_checktype(luaSt, 1, LUA_TUSERDATA);
    return m->func(luaSt, (void*)((char*)lua_touserdata(luaSt, 1) + m->offset));
}

int indexHandler(lua_State* luaSt)
{
    /* stack has userdata, index */
    lua_pushvalue(luaSt, 2);                     /* dup index */
    lua_rawget(luaSt, lua_upvalueindex(1));      /* lookup member by name */
    if (!lua_islightuserdata(luaSt, -1))
    {
        lua_pop(luaSt, 1);                         /* drop value */
        lua_pushvalue(luaSt, 2);                   /* dup index */
        lua_gettable(luaSt, lua_upvalueindex(2));  /* else try methods */
        if (lua_isnil(luaSt, -1))                  /* invalid member */
            luaL_error(luaSt, "cannot get member '%s'", lua_tostring(luaSt, 2));
        return 1;
    }

    return Xet_call(luaSt);                      /* call get function */
}

int newIndexHandler(lua_State* luaSt)
{
    /* stack has userdata, index, value */
    lua_pushvalue(luaSt, 2);                     /* dup index */
    lua_rawget(luaSt, lua_upvalueindex(1));      /* lookup member by name */
    if (!lua_islightuserdata(luaSt, -1))         /* invalid member */
        luaL_error(luaSt, "cannot set member '%s'", lua_tostring(luaSt, 2));

    return Xet_call(luaSt);                      /* call set function */
}

int registerStruct(lua_State* luaSt,
                   const char* name,
                   const Xet_reg_pre* getters,
                   const Xet_reg_pre* setters,
                   const luaL_reg* regMetaMethods,
                   const luaL_reg* regMethods)
{
    int metatable, methods;

    /* create methods table, & add it to the table of globals */
    luaL_openlib(luaSt, name, regMethods, 0);
    methods = lua_gettop(luaSt);

    /* create metatable for NBodyCtx, & add it to the registry */
    luaL_newmetatable(luaSt, name);
    luaL_openlib(luaSt, 0, regMetaMethods, 0);  /* fill metatable */
    metatable = lua_gettop(luaSt);

    lua_pushliteral(luaSt, "__metatable");
    lua_pushvalue(luaSt, methods);    /* dup methods table*/
    lua_rawset(luaSt, metatable);     /* hide metatable:
                                         metatable.__metatable = methods */
    lua_pushliteral(luaSt, "__index");
    lua_pushvalue(luaSt, metatable);     /* upvalue index 1 */
    Xet_add(luaSt, getters);             /* fill metatable with getters */
    lua_pushvalue(luaSt, methods);       /* upvalue index 2 */
    lua_pushcclosure(luaSt, indexHandler, 2);
    lua_rawset(luaSt, metatable);        /* metatable.__index = indexHandler */

    lua_pushliteral(luaSt, "__newindex");
    lua_newtable(luaSt);                 /* table for members you can set */
    Xet_add(luaSt, setters);             /* fill with setters */
    lua_pushcclosure(luaSt, newIndexHandler, 1);
    lua_rawset(luaSt, metatable);     /* metatable.__newindex = newIndexHandler */

    lua_pop(luaSt, 1);            /* drop metatable */
    return 1;                     /* return methods on the stack */
}

int pushEnum(lua_State* luaSt, const MWEnumAssociation* table, int val)
{
    const MWEnumAssociation* p = table;

    while (p->enumVal != -1 && p->enumVal != val)
        ++p;

    if (p->enumVal == -1)
        luaL_error(luaSt, "Got invalid enum value %d", val);

    lua_pushstring(luaSt, p->enumName);

    return 1;
}

int readEnumFromString(lua_State* luaSt, const MWEnumAssociation* table)
{
    const char* str;
    const MWEnumAssociation* p = table;

    str = luaL_checklstring(luaSt, 1, NULL);
    while (p->enumName)
    {
        if (!strcasecmp(p->enumName, str))
            return p->enumVal;
        ++p;
    }

    return InvalidEnum;
}


#if 0
/* FIXME: Source of seed, generator */
static int luaGeneratePlummer(lua_State* luaSt)
{
    dsfmt_t dsfmtState;
    long seed;

    seed = luaL_checklong(luaSt, 1);

    //dsfmt_init_gen_rand(&dsfmtState, seed);
    //generatePlummer(&dsfmtState, );

    return 1;
}
#endif


