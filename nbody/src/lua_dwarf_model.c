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


#include <lua.h>
#include <lauxlib.h>

#include "nbody_types.h"
#include "show.h"
#include "lua_type_marshal.h"
#include "lua_dwarf_model.h"
#include "lua_initial_conditions.h"

#include "milkyway_util.h"

DwarfModel* checkDwarfModel(lua_State* luaSt, int index)
{
    return (DwarfModel*) mw_checknamedudata(luaSt, index, DWARF_MODEL_TYPE);
}

int pushDwarfModel(lua_State* luaSt, const DwarfModel* dm)
{
    DwarfModel* ldm;

    ldm = (DwarfModel*) lua_newuserdata(luaSt, sizeof(DwarfModel));
    if (!ldm)
    {
        warn("Creating DwarfModel userdata failed\n");
        return 1;
    }

    *ldm = *dm;

    //lua_pushvalue(luaSt, -1); /* Copy since luaL_ref() pops */
    //ldm->objRef = luaL_ref(luaSt, LUA_REGISTRYINDEX);
    ldm->generator = 0;

    luaL_getmetatable(luaSt, DWARF_MODEL_TYPE);
    lua_setmetatable(luaSt, -2);

    /* Give this object a new function environment */
    lua_newtable(luaSt);
    lua_setfenv(luaSt, -2);

    return 0;
}

static const MWEnumAssociation dwarfModelOptions[] =
{
    { "Other",      DwarfModelOther   },
    { "Plummer",    DwarfModelPlummer },
    { "King",       DwarfModelKing    },
    { "Dehnen",     DwarfModelDehnen  },
    END_MW_ENUM_ASSOCIATION
};

static int createDwarfModel(lua_State* luaSt)
{
    DwarfModel d = EMPTY_DWARF_MODEL;

    warn("Creating dwarf model\n");

    d.type = checkEnum(luaSt, dwarfModelOptions, -1);
    pushDwarfModel(luaSt, &d);
    return 1;
}

int getDwarfModelT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, dwarfModelOptions, *(int*) v);
}

static int toStringDwarfModel(lua_State* luaSt)
{
    DwarfModel* d;
    char* str;

    d = checkDwarfModel(luaSt, 1);
    str = showDwarfModel(d);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static int gcDwarfModel(lua_State* luaSt)
{
    DwarfModel* dm;

    printf("GOODBYE: Collecting dwarf model\n");
    dm = (DwarfModel*) lua_touserdata(luaSt, 1);
    luaL_unref(luaSt, LUA_REGISTRYINDEX, dm->generator);

    return 0;
}

static const luaL_reg metaMethodsDwarfModel[] =
{
    { "__tostring", toStringDwarfModel },
    { "__gc",       gcDwarfModel       },
    { NULL, NULL }
};

static const luaL_reg methodsDwarfModel[] =
{
    { "create",   createDwarfModel },
    { NULL, NULL }
};

static const Xet_reg_pre gettersDwarfModel[] =
{
    { "type",              getDwarfModelT,       offsetof(DwarfModel, type)              },
    { "nbody" ,            getInt,               offsetof(DwarfModel, nbody)             },
    { "mass" ,             getNumber,            offsetof(DwarfModel, mass)              },
    { "ignoreFinal",       getBool,              offsetof(DwarfModel, ignoreFinal)       },
    { "generator",         getLuaClosure,        offsetof(DwarfModel, generator)         },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersDwarfModel[] =
{
  //{ "type",              setDwarfModelT,       offsetof(DwarfModel, type)              },
    { "nbody" ,            setInt,               offsetof(DwarfModel, nbody)             },
    { "mass" ,             setNumber,            offsetof(DwarfModel, mass)              },
    { "ignoreFinal",       setBool,              offsetof(DwarfModel, ignoreFinal)       },
    { "generator",         setLuaClosure,        offsetof(DwarfModel, generator)         },
    { NULL, NULL, 0 }
};

int registerDwarfModel(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          DWARF_MODEL_TYPE,
                          gettersDwarfModel,
                          settersDwarfModel,
                          metaMethodsDwarfModel,
                          methodsDwarfModel);
}

