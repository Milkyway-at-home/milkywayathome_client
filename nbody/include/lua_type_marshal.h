/*
Copyright (C) 2011  Matthew Arsenault

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _LUA_TYPE_MARSHAL_H_
#define _LUA_TYPE_MARSHAL_H_

#include <lua.h>
#include <lauxlib.h>

typedef int (*Xet_func) (lua_State* luaSt, void* v);

/* member info for get and set handlers */
typedef const struct
{
    const char* name;  /* member name */
    Xet_func func;     /* get or set function for type of member */
    size_t offset;     /* offset of member within NBodyCtx */
}  Xet_reg_pre;

typedef Xet_reg_pre* Xet_reg;

typedef struct
{
    const char* enumName;
    int enumVal;
} MWEnumAssociation;



int getInt(lua_State* luaSt, void* v);
int setInt(lua_State* luaSt, void* v);

int getNumber(lua_State* luaSt, void* v);
int setNumber(lua_State* luaSt, void* v);

int getString(lua_State* luaSt, void* v);


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
int checkEnum(lua_State* luaSt, const MWEnumAssociation* table, int index);

int mw_lua_checkboolean(lua_State* luaSt, int index);
void* mw_checknamedudata(lua_State* luaSt, int index, const char* typeName);

int getBool(lua_State* luaSt, void* v);
int setBool(lua_State* luaSt, void* v);
int getLong(lua_State* luaSt, void* v);
int setLong(lua_State* luaSt, void* v);

#endif /* _LUA_TYPE_MARSHAL_H_ */

