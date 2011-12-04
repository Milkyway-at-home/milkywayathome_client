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

#if !defined(_MILKYWAY_LUA_H_INSIDE_) && !defined(MILKYWAY_LUA_COMPILATION)
  #error "Only milkyway_lua.h can be included directly."
#endif

#ifndef _MILKYWAY_LUA_MARSHAL_H_
#define _MILKYWAY_LUA_MARSHAL_H_

#include <lua.h>
#include <lauxlib.h>

#include "milkyway_extra.h"
#include "milkyway_math.h"

/* All Lua allocations will be aligned to this */
#define MW_LUA_ALIGN 16



typedef int (*Xet_func) (lua_State* luaSt, void* v);


typedef char* (*StructShowFunc) (void*);
typedef const char* (*EnumShowFunc) (int);
typedef void* (*LuaTypeCheckFunc) (lua_State* luaSt, int idx);


/* member info for get and set handlers */
typedef const struct
{
    const char* name;  /* member name */
    Xet_func func;     /* get or set function for type of member */
    size_t offset;     /* offset of member within a struct */
}  Xet_reg_pre;

typedef Xet_reg_pre* Xet_reg;

typedef struct
{
    const char* enumName;
    int enumVal;
} MWEnumAssociation;

#define END_MW_ENUM_ASSOCIATION { NULL, InvalidEnum }

typedef struct
{
    const char* name;
    int type;
    const char* userDataTypeName;
    mwbool required;
    void* value;
} MWNamedArg;

#define END_MW_NAMED_ARG { NULL, -1, NULL, FALSE, NULL }

#define mw_lua_assert_top_type(luaSt, t) assert(lua_type((luaSt), -1) == (t))

int oneTableArgument(lua_State* luaSt, const MWNamedArg* argTable);

int getInt(lua_State* luaSt, void* v);
int setInt(lua_State* luaSt, void* v);

int getLong(lua_State* luaSt, void* v);
int setLong(lua_State* luaSt, void* v);

int getBool(lua_State* luaSt, void* v);
int setBool(lua_State* luaSt, void* v);

int getNumber(lua_State* luaSt, void* v);
int setNumber(lua_State* luaSt, void* v);

int getString(lua_State* luaSt, void* v);

int getCClosure0(lua_State* luaSt, void* v);
int getCClosure1(lua_State* luaSt, void* v);
int getCClosure2(lua_State* luaSt, void* v);

int setCClosure(lua_State* luaSt, void* v);

int setLuaClosure(lua_State* luaSt, void* v);
int getLuaClosure(lua_State* luaSt, void* ref);



void Xet_add(lua_State* luaSt, Xet_reg l);
int Xet_call(lua_State* luaSt);

int indexHandler(lua_State* luaSt);
int newIndexHandler(lua_State* luaSt);

int registerStruct(lua_State* luaSt,
                   const char* name,
                   const Xet_reg_pre* getters,
                   const Xet_reg_pre* setters,
                   const luaL_reg* regMetaMethods,
                   const luaL_reg* regMethods);

int pushEnum(lua_State* luaSt, const MWEnumAssociation* table, int val);
int readEnum(lua_State* luaSt, const MWEnumAssociation* options, const char* name);
int checkEnum(lua_State* luaSt, const MWEnumAssociation* table, int idx);
int expectEnum(lua_State* luaSt, const MWEnumAssociation* table, int idx);


int mw_lua_getglobalfunction(lua_State* luaSt, const char* name);

int mw_lua_checkboolean(lua_State* luaSt, int idx);
mwbool mw_lua_optboolean(lua_State* luaSt, int nArg, mwbool def);

int mw_lua_checktable(lua_State* luaSt, int idx);

lua_CFunction mw_lua_checkcclosure(lua_State* luaSt, int idx);
int mw_lua_checkluaclosure(lua_State* luaSt, int idx);
void mw_lua_pushluaclosure(lua_State* luaSt, int ref);

int mw_lua_checkfunction(lua_State* luaSt, int idx);

void* mw_checknamedudata(lua_State* luaSt, int idx, const char* typeName);
void* mw_tonamedudata(lua_State* luaSt, int ud, const char* typeName);
int mw_lua_typecheck(lua_State* luaSt, int idx, int expectedType, const char* typeName);

void handleNamedArgumentTable(lua_State* luaSt, const MWNamedArg* args, int table);

int pushRealArray(lua_State* luaSt, const real* arr, int n);
real* popRealArray(lua_State* luaSt, int* outN);

int expectTable(lua_State* luaSt, int idx);

int toStringType(lua_State* luaSt, StructShowFunc show, LuaTypeCheckFunc checker);

int pushType(lua_State* luaSt, const char* typeName, size_t typeSize, void* p);
int pushLightType(lua_State* luaSt, const char* typeName, void* p);

void* expectType(lua_State* luaSt, int idx, const char* typeName);

void setModelTableItem(lua_State* luaSt, int table, lua_CFunction generator, const char* name);

#endif /* _MILKYWAY_LUA_MARSHAL_H_ */

