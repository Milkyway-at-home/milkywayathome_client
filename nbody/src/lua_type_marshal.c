
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include <lua.h>
#include <lauxlib.h>

#include "lua_type_marshal.h"
#include "nbody_types.h"
#include "milkyway_util.h"

int mw_lua_checkboolean(lua_State* luaSt, int index)
{
    if (!lua_isboolean(luaSt, index))
        luaL_typerror(luaSt, index, "boolean");

    return lua_toboolean(luaSt, index);
}

lua_CFunction mw_lua_checkcclosure(lua_State* luaSt, int index)
{
    if (!lua_iscfunction(luaSt, index))
        luaL_typerror(luaSt, index, "cclosure");

    return lua_tocfunction(luaSt, index);
}

/* Return reference to Lua closure at index */
int mw_lua_checkluaclosure(lua_State* luaSt, int index)
{
    /* LUA_TFUNCTION can refer to either a C function from the API
     * side, or a Lua closure. lua_tocfunction() returns NULL if it's a Lua function. */

    if (lua_iscfunction(luaSt, index) && !lua_tocfunction(luaSt, index))
        luaL_typerror(luaSt, index, "Lua closure");

    /* Copy since luaL_ref pops and no other lua_check* functions change the stack */
    lua_pushvalue(luaSt, -1);
    return luaL_ref(luaSt, LUA_REGISTRYINDEX);
}

void* mw_checknamedudata(lua_State* luaSt, int index, const char* typeName)
{
    void* v;
    char buf[128];

    if (snprintf(buf, sizeof(buf), "`%s' expected", typeName) == sizeof(buf))
         mw_panic("Error message buffer too small for expected type name\n");

    v = luaL_checkudata(luaSt, index, typeName);
    luaL_argcheck(luaSt, v != NULL, index, buf);

    return v;
}

/* Something between lua_touserdata and lua_checkudata. Checks that it's userdata and
 * of the correct name. Return NULL on failure, rather than erroring. */
void* mw_tonamedudata(lua_State* luaSt, int ud, const char* typeName)
{
    void* p;

    p = lua_touserdata(luaSt, ud);
    if (!p)
        return NULL;

    lua_getfield(luaSt, LUA_REGISTRYINDEX, typeName);  /* get correct metatable */
    if (!lua_getmetatable(luaSt, ud))
    {
        lua_pop(luaSt, 1);
        return NULL;
    }

    if (!lua_rawequal(luaSt, -1, -2))
        p = NULL;

    lua_pop(luaSt, 2);
    return p;
}

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

int getLong(lua_State* luaSt, void* v)
{
    lua_pushinteger(luaSt, *(long*)v);
    return 1;
}

int setLong(lua_State* luaSt, void* v)
{
    *(long*)v = luaL_checklong(luaSt, 3);
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

int getBool(lua_State* luaSt, void* v)
{
    lua_pushboolean(luaSt, *(int*)v);
    return 1;
}

int setBool(lua_State* luaSt, void* v)
{
    *(int*)v = mw_lua_checkboolean(luaSt, 3);
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

    lua_pop(luaSt, 2);         /* drop metatable and methods */
    return 0;
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

static int checkEnumError(lua_State* luaSt, const MWEnumAssociation* p, const char* badStr)
{
    const MWEnumAssociation* nextP;
    char errBuf[2048] = "Expected enum value where options are: ";
    char badOpt[1024] = "";
    size_t badSize, enumLen, errLen;

    errLen = strlen(errBuf);
    while (p->enumName)
    {
        nextP = &p[1];
        enumLen = strlen(p->enumName);
        if (errLen + enumLen + 6 > sizeof(errBuf))  /* max possible */
            mw_panic("Enum options too large for error string buffer!\n");

        errLen += enumLen;
        if (nextP->enumName) /* If there is a next one, use a comma */
        {
            strcat(strcat(errBuf, p->enumName), ", ");
            errLen += 2;  /* ', ' */
        }
        else
        {
            strcat(strcat(strcat(errBuf, "or "), p->enumName), ".");
            errLen += 4; /* 'or ' + '.' */
        }

        p = nextP;
    }

    /* If there's extra space, might as well say what the bad option was */
    badSize = snprintf(badOpt, sizeof(badOpt), " Invalid option '%s'", badStr);
    if (   (badSize != sizeof(badOpt))
        && (badSize < (sizeof(errBuf) - errLen + 3)))
    {
        strncat(errBuf, badOpt, sizeof(errBuf));
    }

    return luaL_argerror(luaSt, 1, errBuf);
}

int checkEnum(lua_State* luaSt, const MWEnumAssociation* table, int index)
{
    const char* str;
    const MWEnumAssociation* p = table;

    str = luaL_checklstring(luaSt, index, NULL);

    while (p->enumName)
    {
        if (!strcasecmp(p->enumName, str))
            break;
        ++p;
    }

    if (!p->enumName)
    {
        checkEnumError(luaSt, table, str);
        return InvalidEnum;
    }

    return p->enumVal;
}


static void setValueFromType(lua_State* luaSt, void* v, int type, int index)
{
    switch (type)
    {
        case LUA_TNUMBER:
            *(real*) v = (real) lua_tonumber(luaSt, index);
            break;

        case LUA_TBOOLEAN:
            *(mwbool*) v = (mwbool) lua_toboolean(luaSt, index);
            break;

        case LUA_TSTRING:
            *(char**) v = strdup(lua_tostring(luaSt, index));
            break;

        case LUA_TUSERDATA:
            *(void**) v = lua_touserdata(luaSt, index);
            break;

        case LUA_TTABLE:
        case LUA_TFUNCTION:

        case LUA_TLIGHTUSERDATA:

        case LUA_TTHREAD:
        case LUA_TNIL:
        default:
            mw_panic("Unhandled type %s (%d)\n", luaL_typename(luaSt, type), type);
    }
}

static void namedArgumentError(lua_State* luaSt, const MWNamedArg* p, int arg, int index)
{
    luaL_error(luaSt, "Bad argument for key '%s' in argument #%d (`%s' expected, got %s)",
               p->name,
               arg,
               p->userDataTypeName ? p->userDataTypeName : lua_typename(luaSt, p->luaType),
               lua_type(luaSt, index) == LUA_TUSERDATA ? "other userdata" : luaL_typename(luaSt, index)
        );
}

/* Check if any keys are left in the table and error if there are */
static void checkExtraArguments(lua_State* luaSt, int table)
{
    unsigned int extraCount = 0;

    lua_pushnil(luaSt);  /* first key */
    while (lua_next(luaSt, table) != 0) /* Iterate remaining keys */
    {
        ++extraCount;
        lua_pop(luaSt, 1);
    }

    /* TODO: Print bad keys, better message */
    if (extraCount > 0)
        luaL_error(luaSt, "%d unknown arguments found", extraCount);
}

void handleNamedArgumentTable(lua_State* luaSt, const MWNamedArg* args, int table)
{
    const MWNamedArg* p;
    char buf[128];
    int type, item;

    p = args;
    while (p->name)
    {
        lua_pushstring(luaSt, p->name);
        lua_pushvalue(luaSt, -1);  /* Copy the key since lua_gettable pops it */

        lua_gettable(luaSt, table);
        item = lua_gettop(luaSt);

        if (lua_isnil(luaSt, item) && p->required)
        {
            if (snprintf(buf, sizeof(buf), "Missing required named argument '%s'", p->name) == sizeof(buf))
                mw_panic("Error message buffer too small for key name '%s'\n", p->name);

            luaL_argerror(luaSt, 1, buf);
        }


        /* We do our own type checking and errors to avoid
           Confusing and innaccurate error messages, which suggest the use of the table is wrong. */
        type = lua_type(luaSt, -1);
        if (type != p->luaType)
        {
            namedArgumentError(luaSt, p, 1, item);
        }

        if (type == LUA_TUSERDATA) /* We must do another level of checking for the actual type */
        {
            if (!mw_tonamedudata(luaSt, item, p->userDataTypeName))
                namedArgumentError(luaSt, p, 1, -1);
        }
        else if (type == LUA_TLIGHTUSERDATA)
        {
            mw_panic("Unhandled named argument type lightuserdata\n");
        }

        setValueFromType(luaSt, p->value, p->luaType, item);

        lua_pop(luaSt, 1);

        /* Top item, should now be at copy of key */
        lua_pushnil(luaSt);
        lua_rawset(luaSt, table); /* Wipe out this element so we can check for unknown arguments */

        ++p;
    }

    checkExtraArguments(luaSt, table);
}

static inline int getCClosureN(lua_State* luaSt, void* v, int n)
{
    lua_pushcclosure(luaSt, (lua_CFunction) v, n);
    return 1;
}

int getCClosure0(lua_State* luaSt, void* v)
{
    return getCClosureN(luaSt, v, 0);
}

int getCClosure1(lua_State* luaSt, void* v)
{
    return getCClosureN(luaSt, v, 1);
}

int getCClosure2(lua_State* luaSt, void* v)
{
    return getCClosureN(luaSt, v, 2);
}

int setCClosure(lua_State* luaSt, void* v)
{
    warn("Setting closure\n");
    *(lua_CFunction*) v = mw_lua_checkcclosure(luaSt, 3);
    return 0;
}

int setLuaClosure(lua_State* luaSt, void* _ref)
{
    int* ref = (int*) _ref;

    if (*ref != LUA_REFNIL && *ref != LUA_NOREF)
        luaL_unref(luaSt, LUA_REGISTRYINDEX, *ref);

    *ref = mw_lua_checkluaclosure(luaSt, 3);
    assert(*ref != LUA_NOREF);

    return 0;
}

int getLuaClosure(lua_State* luaSt, void* ref)
{
    warn("Getting lua closure\n");
    lua_rawgeti(luaSt, LUA_REGISTRYINDEX, *(int*) ref);

    return 1;
}

