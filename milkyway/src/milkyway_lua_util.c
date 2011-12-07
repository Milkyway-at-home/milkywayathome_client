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
#include <lualib.h>
#include <sys/stat.h>

#include "milkyway_lua_marshal.h"
#include "milkyway_util.h"
#include "milkyway_lua_util.h"

#if HAVE_DIRECT_H
  #include <direct.h>
#endif

#if HAVE_SYS_STAT_H
  #include <sys/stat.h>
#endif


#ifdef _WIN32
  /* WTF? Why is the windows one missing the 2nd argument? */
  #define mkdir(x, y) mkdir(x)
#endif

/* Lazy binding for creating directories. Should probably error and things */
static int lua_mkdir(lua_State* luaSt)
{
    const char* path;

    path = luaL_checkstring(luaSt, 1);
    lua_pushinteger(luaSt, mkdir(path, 0777));

    return 1;
}

/* Map over values in a table */
static int luaMap(lua_State* luaSt)
{
    int i, n, f = 1, table = 2, newTable = 3;

    if (lua_gettop(luaSt) != 2)
    {
        luaL_argerror(luaSt, 2, "Expected 2 arguments");
    }

    mw_lua_checktable(luaSt, table);
    mw_lua_checkfunction(luaSt, f);

    n = luaL_getn(luaSt, table);

    lua_newtable(luaSt);
    for (i = 0; i < n; ++i)
    {
        lua_pushvalue(luaSt, f);
        lua_rawgeti(luaSt, table, i + 1);
        lua_call(luaSt, 1, 1);
        lua_rawseti(luaSt, newTable, i + 1);
    }

    return 1;
}

/* FIXME?: can be side effecting on accumulator */
static int luaFoldl(lua_State* luaSt)
{
    int i, n, acc, f = 1, ini = 2, table = 3;

    if (lua_gettop(luaSt) != 3)
        luaL_argerror(luaSt, 3, "Expected 3 arguments");

    mw_lua_checkfunction(luaSt, f);
    /* accumulator can be anything */
    mw_lua_checktable(luaSt, table);

    n = luaL_getn(luaSt, table);

    lua_pushvalue(luaSt, ini);
    acc = lua_gettop(luaSt);

    for (i = 0; i < n; ++i)
    {
        lua_pushvalue(luaSt, f);
        lua_pushvalue(luaSt, acc);
        lua_rawgeti(luaSt, table, i + 1);
        lua_call(luaSt, 2, 1);
        lua_replace(luaSt, acc);
    }

    return 1;
}

static int luaZipWith(lua_State* luaSt)
{
    int i, n, f = 1, tableA = 2, tableB = 3, newTable = 4;

    if (lua_gettop(luaSt) != 3)
        luaL_argerror(luaSt, 3, "Expected 3 arguments");

    mw_lua_checkfunction(luaSt, f);
    mw_lua_checktable(luaSt, tableA);
    mw_lua_checktable(luaSt, tableB);

    n = MIN(luaL_getn(luaSt, tableA), luaL_getn(luaSt, tableB));

    lua_newtable(luaSt);
    for (i = 0; i < n; ++i)
    {
        lua_pushvalue(luaSt, f);
        lua_rawgeti(luaSt, tableA, i + 1);
        lua_rawgeti(luaSt, tableB, i + 1);
        lua_call(luaSt, 2, 1);
        lua_rawseti(luaSt, newTable, i + 1);
    }

    return 1;
}

static int luaGetTime(lua_State* luaSt)
{
    lua_Number t = mwGetTime();
    lua_pushnumber(luaSt, t);
    return 1;
}

int registerUtilityFunctions(lua_State* luaSt)
{
    lua_register(luaSt, "map", luaMap);
    lua_register(luaSt, "foldl", luaFoldl);
    lua_register(luaSt, "zipWith", luaZipWith);
    lua_register(luaSt, "getTime", luaGetTime);
    lua_register(luaSt, "mkdir", lua_mkdir);

    return 0;
}

static const luaL_Reg normalLibs[] =
{
    { "",              luaopen_base   },
    { LUA_TABLIBNAME,  luaopen_table  },
    { LUA_STRLIBNAME,  luaopen_string },
 // { LUA_MATHLIBNAME, luaopen_math   },  We replace math with bindings to whatever math we're using
    { NULL, NULL}
};

static const luaL_Reg debugOnlyLibs[] =
{
    { LUA_LOADLIBNAME, luaopen_package },
    { LUA_IOLIBNAME,   luaopen_io      },
    { LUA_OSLIBNAME,   luaopen_os      },
    { LUA_DBLIBNAME,   luaopen_debug   },
    { LUA_MATHLIBNAME, luaopen_math    },
    { NULL, NULL}
};

static void mw_luaL_register(lua_State* luaSt, const luaL_Reg* libs)
{
    while (libs->name)
    {
        lua_pushcfunction(luaSt, libs->func);
        lua_pushstring(luaSt, libs->name);
        lua_call(luaSt, 1, 0);
        ++libs;
    }
}

/* Finer control over which standard libraries are opened for sandboxing */
void mw_lua_openlibs(lua_State* luaSt, mwbool debug)
{
    mw_luaL_register(luaSt, normalLibs);
    if (debug)
        mw_luaL_register(luaSt, debugOnlyLibs);
}

static int doWithArgs(lua_State* luaSt, const char** args, unsigned int nArgs)
{
    unsigned int i;

    if (args)
    {
        for (i = 0; i < nArgs; ++i)
            lua_pushstring(luaSt, args[i]);
    }

    if (lua_pcall(luaSt, nArgs, 0, 0))
    {
        return 1;
    }

    return 0;
}

int dostringWithArgs(lua_State* luaSt,
                     const char* str,
                     const char** args,
                     unsigned int nArgs)
{
    /* If either fails there will be 1 error on the stack */
    return luaL_loadstring(luaSt, str) || doWithArgs(luaSt, args, nArgs);
}

int dofileWithArgs(lua_State* luaSt,
                   const char* filename,
                   const char** args,
                   unsigned int nArgs)
{
    return luaL_loadfile(luaSt, filename) || doWithArgs(luaSt, args, nArgs);
}

int mwBindBOINCStatus(lua_State* luaSt)
{
    lua_pushboolean(luaSt, BOINC_APPLICATION);
    lua_setglobal(luaSt, "isBOINCApplication");

    lua_pushboolean(luaSt, mw_is_standalone());
    lua_setglobal(luaSt, "isStandalone");

    return 0;
}

/*
  Derived from _mingw_aligned_malloc and friends, implemented using Microsoft's public
  interfaces and with the help of the algorithm description provided
  by Wu Yongwei: http://sourceforge.net/mailarchive/message.php?msg_id=3847075

  I hereby place this implementation in the public domain.
  -- Steven G. Johnson (stevenj@alum.mit.edu)
*/
#define POWER_OF_TWO(n) (!((n) & ((n) - 1)))
#define UI(p) ((uintptr_t) (p))
#define CP(p) ((char*) p)

#define PTR_ALIGN(p0, alignment)                \
    ((void*) (((UI(p0) + (alignment + sizeof(void*)))   \
                & (~UI(alignment - 1)))))

/* Pointer must sometimes be aligned; assume sizeof(void*) is a power of two. */
#define ORIG_PTR(p) (*(((void**) (UI(p) & (~UI(sizeof(void*) - 1)))) - 1))


static void* mw_lua_aligned_malloc(size_t size)
{
    void* p;
    void* p0;

    assert(POWER_OF_TWO(MW_LUA_ALIGN));
    assert(MW_LUA_ALIGN >= sizeof(void*));

    /* Including the extra sizeof(void*) is overkill on a 32-bit
       machine, since malloc is already 8-byte aligned, as long
       as we enforce alignment >= 8 ...but oh well.  */

    p0 = malloc(size + MW_LUA_ALIGN + sizeof(void*));
    if (!p0)
    {
        return NULL;
    }

    p = PTR_ALIGN(p0, MW_LUA_ALIGN);
    ORIG_PTR(p) = p0;
    return p;
}

static void mw_lua_aligned_free(void* memblock)
{
    if (memblock)
    {
        free(ORIG_PTR(memblock));
    }
}

static void* mw_lua_aligned_realloc(void* memblock, size_t size)
{
    void* p;
    void* p0;
    ptrdiff_t shift;

    if (!memblock)
    {
        return mw_lua_aligned_malloc(size);
    }

    if (size == 0)
    {
        mw_lua_aligned_free(memblock);
        return NULL;
    }

    p0 = ORIG_PTR(memblock);

    assert(memblock == PTR_ALIGN(p0, MW_LUA_ALIGN));    /* It is an error for the alignment to change. */
    shift = CP(memblock) - CP(p0);

    p0 = realloc(p0, size + MW_LUA_ALIGN + sizeof(void*));
    if (!p0)
    {
        return NULL;
    }

    p = PTR_ALIGN(p0, MW_LUA_ALIGN);

    /* Relative shift of actual data may be different from before, ugh.  */
    if (shift != CP(p) - CP(p0))
    {
        /* ugh, moves more than necessary if size is increased.  */
        memmove(CP(p), CP(p0) + shift, size);
    }

    ORIG_PTR(p) = p0;
    return p;
}

static void* mwLuaAlignedAllocator(void* ud, void* ptr, size_t osize, size_t nsize)
{
    (void) ud, (void) osize;

    if (nsize == 0)
    {
        mw_lua_aligned_free(ptr);
        return NULL;
    }

    else
    {
        return mw_lua_aligned_realloc(ptr, nsize);
    }

}


/* Add the standard panic function back to our state with the custom allocator */
static int mw_lua_panic(lua_State* luaSt)
{
    (void) luaSt;
    fprintf(stderr,
            "PANIC: unprotected error in call to Lua API (%s)\n",
            lua_tostring(luaSt, -1));
    return 0;
}

lua_State* mw_lua_newstate(void)
{
    lua_State* luaSt;

    luaSt = lua_newstate(mwLuaAlignedAllocator, NULL);
    if (luaSt)
    {
        lua_atpanic(luaSt, mw_lua_panic);
    }

    return luaSt;
}

