/*
 *  Copyright (c) 2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <lua.h>
#include <lauxlib.h>

#include "milkyway_lua_marshal.h"
#include "milkyway_util.h"

/*
  Check if global function exists.
   If symbol exists and is a function, put it on the stack and return 0.
   Return > 0 if wrong type
   Return < 0 if does not exist
*/
int mw_lua_getglobalfunction(lua_State* luaSt, const char* name)
{
    lua_getglobal(luaSt, name);
    if (lua_isnil(luaSt, -1))
    {
        mw_printf("Didn't find required global function '%s'\n", name);
        return -1;
    }
    else if (!lua_isfunction(luaSt, -1))
    {
        mw_printf("Expected required global '%s' to be %s, got %s\n",
                  name, lua_typename(luaSt, LUA_TFUNCTION), luaL_typename(luaSt, -1));
        return 1;
    }
    else
    {
        return 0;
    }
}

int mw_lua_checkboolean(lua_State* luaSt, int idx)
{
    if (!lua_isboolean(luaSt, idx))
        luaL_typerror(luaSt, idx, "boolean");

    return lua_toboolean(luaSt, idx);
}

mwbool mw_lua_optboolean(lua_State* luaSt, int nArg, mwbool def)
{
    return lua_isnoneornil(luaSt, nArg) ? def : mw_lua_checkboolean(luaSt, nArg);
}

/* Check that argument at index is a table. Fails if not, returns the index if it is */
int mw_lua_checktable(lua_State* luaSt, int idx)
{
    return lua_istable(luaSt, idx) ? idx : luaL_typerror(luaSt, idx, "table");
}

lua_CFunction mw_lua_checkcclosure(lua_State* luaSt, int idx)
{
    if (!lua_iscfunction(luaSt, idx))
        luaL_typerror(luaSt, idx, "cclosure");

    return lua_tocfunction(luaSt, idx);
}

int mw_lua_checkfunction(lua_State* luaSt, int idx)
{
    if (!lua_isfunction(luaSt, idx))
        return luaL_typerror(luaSt, idx, "function");

    return 0;
}

/* Return reference to Lua closure at idx */
int mw_lua_checkluaclosure(lua_State* luaSt, int idx)
{
    /* LUA_TFUNCTION can refer to a cclosure or a Lua closure */
    if (!lua_isfunction(luaSt, idx) || lua_iscfunction(luaSt, idx))
    {
        luaL_typerror(luaSt, idx, "Lua closure");
    }

    /* Copy since luaL_ref pops and no other lua_check* functions change the stack */
    lua_pushvalue(luaSt, -1);
    return luaL_ref(luaSt, LUA_REGISTRYINDEX);
}

void mw_lua_pushluaclosure(lua_State* luaSt, int ref)
{
    lua_pushinteger(luaSt, ref);
    lua_rawget(luaSt, LUA_REGISTRYINDEX);
}

void* mw_checknamedudata(lua_State* luaSt, int idx, const char* typeName)
{
    void* v;
    char buf[128];

    if (snprintf(buf, sizeof(buf), "`%s' expected", typeName) == sizeof(buf))
         mw_panic("Error message buffer too small for expected type name\n");

    v = luaL_checkudata(luaSt, idx, typeName);
    luaL_argcheck(luaSt, v != NULL, idx, buf);

    return v;
}

static int mw_lua_equal_userdata_name(lua_State* luaSt, int idx, const char* typeName)
{
    int equal;

    lua_getfield(luaSt, LUA_REGISTRYINDEX, typeName);  /* get correct metatable */
    if (!lua_getmetatable(luaSt, idx))
    {
        mw_printf("Failed to get metatable from index %d\n", idx);
        lua_pop(luaSt, 1);
        return 0;
    }

    equal = lua_rawequal(luaSt, -1, -2);
    lua_pop(luaSt, 2);

    return equal;
}

/* Something between lua_touserdata and lua_checkudata. Checks that it's userdata and
 * of the correct name. Return NULL on failure, rather than erroring. */
void* mw_tonamedudata(lua_State* luaSt, int ud, const char* typeName)
{
    void* p;

    p = lua_touserdata(luaSt, ud);
    if (!p)
        return NULL;

    return mw_lua_equal_userdata_name(luaSt, ud, typeName) ? p : NULL;
}

/* The luaL_check* family of functions are intended for use for Lua
 * stuff written in the C api, not really getting results from
 * lua. luaL_error() and co. abort the whole program, which are caught
 * when running things with lua_pcall, but not what we want for
 * getting what we want with within C side stuff.

 This is for printing an error when manually checking values returned
 to C stuff. Returns 1 if type problem, 0 otherwise and places appropriate error on stack */
int mw_lua_typecheck(lua_State* luaSt, int idx, int expectedType, const char* typeName)
{
    int type;

    type = lua_type(luaSt, idx);
    if (   type == expectedType  /* got userdata of wrong type */
        && expectedType == LUA_TUSERDATA
        && !mw_lua_equal_userdata_name(luaSt, idx, typeName))
    {
        /* TODO: Get typename of wong userdata type, which is kind of the point of checking for this */
        lua_pushfstring(luaSt, "Type error: userdata %s expected, got other userdata", typeName);
        return 1;
    }
    else if (type != expectedType) /* Anything else is wrong. */
    {
        lua_pushfstring(luaSt,
                        "Type error: %s expected, got %s",
                        lua_typename(luaSt, expectedType),
                        luaL_typename(luaSt, idx));
        return 1;
    }
    else
    {
        return 0;
    }
}

int oneTableArgument(lua_State* luaSt, const MWNamedArg* argTable)
{
    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 table argument");

    handleNamedArgumentTable(luaSt, argTable, 1);
    return 0;
}

/* Mostly from example at http://lua-users.org/wiki/BindingWithMembersAndMethods */


int getInt(lua_State* luaSt, void* v)
{
    lua_pushnumber(luaSt, *(int*) v);
    return 1;
}

int setInt(lua_State* luaSt, void* v)
{
    *(int*) v = luaL_checkint(luaSt, 3);
    return 0;
}

int getLong(lua_State* luaSt, void* v)
{
    lua_pushinteger(luaSt, *(long*) v);
    return 1;
}

int setLong(lua_State* luaSt, void* v)
{
    *(long*) v = luaL_checklong(luaSt, 3);
    return 0;
}

int getNumber(lua_State* luaSt, void* v)
{
    lua_pushnumber(luaSt, (lua_Number) *(real*) v);
    return 1;
}

int setNumber(lua_State* luaSt, void* v)
{
    *(real*) v = (real) luaL_checknumber(luaSt, 3);
    return 0;
}

int getBool(lua_State* luaSt, void* v)
{
    lua_pushboolean(luaSt, *(int*) v);
    return 1;
}

int setBool(lua_State* luaSt, void* v)
{
    *(int*) v = mw_lua_checkboolean(luaSt, 3);
    return 0;
}

int getString(lua_State* luaSt, void* v)
{
    lua_pushstring(luaSt, (char*) v);
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

    /* create methods table, and add it to the table of globals */
    luaL_register(luaSt, name, regMethods);
    methods = lua_gettop(luaSt);

    /* create metatable for type, and add it to the registry */
    luaL_newmetatable(luaSt, name);
    luaL_register(luaSt, NULL, regMetaMethods);  /* fill metatable */
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

static void checkEnumErrorStr(char* errBuf, size_t errBufSize, const MWEnumAssociation* p, const char* badStr)
{
    const MWEnumAssociation* nextP;
    static const char errStart[] = "Expected enum value where options are: ";
    char badOpt[1024];
    size_t badSize, enumLen, errLen;
    size_t remSize; /* Remaining size in buffer */

    strcpy(errBuf, errStart);
    errLen = sizeof(errStart) - 1;
    while (p->enumName)
    {
        nextP = &p[1];
        enumLen = strlen(p->enumName);
        if (errLen + enumLen + 6 > errBufSize)  /* max possible */
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
    remSize = sizeof(errBuf) - errLen - 1;
    if ((badSize != sizeof(badOpt)) && (badSize < remSize))
    {
        strncat(errBuf, badOpt, remSize);
    }
}

static int findEnumEntry(const MWEnumAssociation* table, const char* str)
{
    const MWEnumAssociation* p = table;

    if (!str)
    {
        return InvalidEnum;
    }

    while (p->enumName)
    {
        if (!strcasecmp(p->enumName, str))
            break;
        ++p;
    }

    return p->enumVal;
}

/* Check for the expected enum value at idx. Put error string on stack */
int expectEnum(lua_State* luaSt, const MWEnumAssociation* table, int idx)
{
    int entry;
    const char* checkStr;
    char errBuf[4096];

    checkStr = lua_tostring(luaSt, idx);
    entry = findEnumEntry(table, checkStr);
    if (entry != InvalidEnum)
    {
        return entry;
    }

    checkEnumErrorStr(errBuf, sizeof(errBuf), table, checkStr);
    mw_printf("%s\n", errBuf);
    lua_pop(luaSt, 1);
    return InvalidEnum;
}

/* Check for the expected enum value at idx. Error if invalid value */
int checkEnum(lua_State* luaSt, const MWEnumAssociation* table, int idx)
{
    int entry;
    const char* str;
    char errBuf[4096];

    str = luaL_checklstring(luaSt, idx, NULL);
    entry = findEnumEntry(table, str);
    if (entry != InvalidEnum)
    {
        return entry;
    }

    checkEnumErrorStr(errBuf, sizeof(errBuf), table, str);
    luaL_argerror(luaSt, 1, errBuf);
    return InvalidEnum;
}


int readEnum(lua_State* luaSt, const MWEnumAssociation* options, const char* name)
{
    int rc;

    assert(name);
    lua_pushstring(luaSt, name);
    rc = checkEnum(luaSt, options, lua_gettop(luaSt));
    lua_pop(luaSt, 1);
    return rc;
}

static void setValueFromType(lua_State* luaSt, void* v, int type, int idx)
{
    switch (type)
    {
        case LUA_TNUMBER:
            *(real*) v = (real) lua_tonumber(luaSt, idx);
            break;

        case LUA_TBOOLEAN:
            *(mwbool*) v = (mwbool) lua_toboolean(luaSt, idx);
            break;

        case LUA_TSTRING:
            *(const char**) v = lua_tostring(luaSt, idx);
            break;

        case LUA_TUSERDATA:
            *(void**) v = lua_touserdata(luaSt, idx);
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

static void namedArgumentError(lua_State* luaSt, const MWNamedArg* p, int arg, int idx)
{
    luaL_error(luaSt, "Bad argument for key '%s' in argument #%d (`%s' expected, got %s)",
               p->name,
               arg,
               p->userDataTypeName ? p->userDataTypeName : lua_typename(luaSt, p->type),
               lua_type(luaSt, idx) == LUA_TUSERDATA ? "other userdata" : luaL_typename(luaSt, idx)
        );
}

/* Lazy search of the names table */
static mwbool nameInArgTable(const MWNamedArg* p, const char* key)
{
    if (!key)
        return FALSE;

    while (p->name)
    {
        if (!strcmp(p->name, key))
            return TRUE;
        ++p;
    }

    return FALSE;
}

static int pushUnknownNamedArgumentError(lua_State* luaSt, int keyIdx, mwbool typeError)
{
    if (typeError)
    {
        lua_pushfstring(luaSt,
                        "  Unexpected named argument of non-string type '%s'\n",
                        luaL_typename(luaSt, keyIdx));
    }
    else
    {
        lua_pushfstring(luaSt, "  Unknown named argument '%s'\n", lua_tostring(luaSt, keyIdx));
    }

    return 1;
}

static int pushMultilineErrorPrefix(lua_State* luaSt)
{
    luaL_where(luaSt, 1);
    lua_pushliteral(luaSt, "\n");
    lua_concat(luaSt, 2);

    return 1;
}

/* Check if any keys in the table aren't in the expected table, and error if there are */
static void checkExtraArguments(lua_State* luaSt, const MWNamedArg* args, int table)
{
    int iniTop, keyIdx, extraCount;
    const char* key = NULL;
    mwbool badFound = FALSE, keyIsString = FALSE, keyIsInTable = FALSE;

    iniTop = lua_gettop(luaSt);

    lua_pushnil(luaSt);  /* first key */
    while (lua_next(luaSt, table) != 0) /* Iterate keys in table */
    {
        lua_pop(luaSt, 1);    /* Get rid of value, keep key on top  */
        keyIdx = lua_gettop(luaSt);

        /* Named arguments can only be strings */
        keyIsString = lua_isstring(luaSt, keyIdx) && !lua_isnumber(luaSt, keyIdx);

        /* Avoid automatic conversions from numbers to string */
        key = keyIsString ? lua_tostring(luaSt, keyIdx) : NULL;

        /* Only look for the key if it's actually a string */
        keyIsInTable = nameInArgTable(args, key);

        if (!keyIsString || !keyIsInTable)  /* Bad argument */
        {
            if (!badFound)  /* Add prefix of filename, line number for first error */
            {
                badFound = TRUE;
                pushMultilineErrorPrefix(luaSt);

                lua_insert(luaSt, iniTop + 1); /* Save this at the bottom of the extra errors for later */
                ++keyIdx; /* The insert just shifted this up 1 place */
            }

            /* Generate appropriate error message depending on whether it's actually a string */
            pushUnknownNamedArgumentError(luaSt, keyIdx, !keyIsString);

            /* Make sure key is still on top for iterating over the table */
            lua_insert(luaSt, keyIdx);
        }
    }

    extraCount = lua_gettop(luaSt) - iniTop - 1; /* One line for each argument + prefix */
    if (extraCount > 0)
    {
        lua_pushfstring(luaSt, "  %d bad named arguments found", extraCount);
        lua_concat(luaSt, extraCount + 2); /* + 1 for total argument count line, + 1 for prefix line */
        lua_error(luaSt);
    }

    assert(iniTop == lua_gettop(luaSt));
}

/* Some of the automatic type conversions are acceptable, others are
 * not.
 * Allowing numbers as strings to work lets us simply give the
 * extra arguments to the script, so we can have non-number
 * arguments.*/
static int typeEqualOrConversionOK(lua_State* luaSt, int type, int idx)
{
    int luaType = lua_type(luaSt, idx);

    if (type == luaType)
        return TRUE;

    if (type == LUA_TNUMBER && lua_isnumber(luaSt, idx))
        return TRUE;

    return FALSE;
}

void handleNamedArgumentTable(lua_State* luaSt, const MWNamedArg* args, int table)
{
    const MWNamedArg* p;
    char buf[128];
    int item;

    p = args;
    while (p->name)
    {
        lua_pushstring(luaSt, p->name);
        lua_pushvalue(luaSt, -1);  /* Copy the key since lua_gettable pops it */

        lua_gettable(luaSt, table);
        item = lua_gettop(luaSt);

        if (lua_isnil(luaSt, item))
        {
            if (!p->required)
            {
                lua_pop(luaSt, 2);
                ++p;
                continue;
            }

            if (snprintf(buf, sizeof(buf), "Missing required named argument '%s'", p->name) == sizeof(buf))
                mw_panic("Error message buffer too small for key name '%s'\n", p->name);

            luaL_argerror(luaSt, table, buf);
        }

        /* We do our own type checking and errors to avoid
           Confusing and innaccurate error messages, which suggest the use of the table is wrong. */
        if (!typeEqualOrConversionOK(luaSt, p->type, -1))
        {
            namedArgumentError(luaSt, p, table, item);
        }

        if (p->type == LUA_TUSERDATA) /* We must do another level of checking for the actual type */
        {
            if (!mw_tonamedudata(luaSt, item, p->userDataTypeName))
                namedArgumentError(luaSt, p, table, -1);
        }
        else if (p->type == LUA_TLIGHTUSERDATA)
        {
            mw_panic("Unhandled named argument type lightuserdata\n");
        }

        setValueFromType(luaSt, p->value, p->type, item);
        lua_pop(luaSt, 2);
        ++p;
    }

    checkExtraArguments(luaSt, args, table);
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
    lua_rawgeti(luaSt, LUA_REGISTRYINDEX, *(int*) ref);

    return 1;
}

int pushRealArray(lua_State* luaSt, const real* arr, int n)
{
    int i, table;

    lua_createtable(luaSt, n, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < n; ++i)
    {
        lua_pushnumber(luaSt, (lua_Number) arr[i]);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 0;
}

real* popRealArray(lua_State* luaSt, int* outN)
{
    real* arr;
    int i, n, table;

    table = lua_gettop(luaSt);
    luaL_checktype(luaSt, table, LUA_TTABLE);
    n = luaL_getn(luaSt, table);  /* get size of table */

    arr = (real*) mwMalloc(sizeof(real) * n);
    for (i = 0; i < n; ++i)
    {
        lua_rawgeti(luaSt, table, i + 1);  /* push t[i] */
        luaL_checktype(luaSt, -1, LUA_TNUMBER);
        arr[i] = lua_tonumber(luaSt, -1);
        lua_pop(luaSt, 1);
    }

    lua_pop(luaSt, 1);

    if (outN)
        *outN = n;

    return arr;
}

int expectTable(lua_State* luaSt, int idx)
{
    return mw_lua_typecheck(luaSt, idx, LUA_TTABLE, NULL);
}

/* Create a function for __tostring metamethods from show* functions */
int toStringType(lua_State* luaSt, StructShowFunc show, LuaTypeCheckFunc checker)
{
    void* p;
    char* str;

    p = checker(luaSt, 1);
    assert(p);
    str = show(p);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

int pushType(lua_State* luaSt, const char* typeName, size_t typeSize, void* p)
{
    void* lp;

    lp = lua_newuserdata(luaSt, typeSize);
    if (!lp)
    {
        mw_printf("Creating userdata '%s' failed\n", typeName);
        return 0;
    }

    assert((uintptr_t) lp % MW_LUA_ALIGN == 0); /* This must be true for dSFMT intrinsics stuff to work */

    luaL_getmetatable(luaSt, typeName);
    lua_setmetatable(luaSt, -2);

#if 0
    /* Give this object a new function environment; for installing
     * arbitrary lua functions into a type */
    lua_newtable(luaSt);
    lua_setfenv(luaSt, -2);
#endif

    memcpy(lp, p, typeSize);

    return 1;
}

/* toType* functions with a bonus error message */
void* expectType(lua_State* luaSt, int idx, const char* typeName)
{
    void* p;

    p = mw_tonamedudata(luaSt, idx, typeName);
    if (!p)
        mw_lua_typecheck(luaSt, idx, LUA_TUSERDATA, typeName);

    return p;
}

void setModelTableItem(lua_State* luaSt, int table, lua_CFunction generator, const char* name)
{
    lua_pushcfunction(luaSt, generator);
    lua_setfield(luaSt, table, name);
}

