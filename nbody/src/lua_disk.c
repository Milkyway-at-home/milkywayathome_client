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
#include "lua_disk.h"

#include "milkyway_util.h"

Disk* checkDisk(lua_State* luaSt, int index)
{
    return (Disk*) mw_checknamedudata(luaSt, index, DISK_TYPE);
}

int pushDisk(lua_State* luaSt, const Disk* d)
{
    Disk* ld;

    ld = (Disk*) lua_newuserdata(luaSt, sizeof(Disk));
    if (!ld)
    {
        warn("Creating Disk userdata failed\n");
        return 1;
    }

    luaL_getmetatable(luaSt, DISK_TYPE);
    lua_setmetatable(luaSt, -2);

    *ld = *d;

    return 0;
}

static const MWEnumAssociation diskOptions[] =
{
    { "exponential",    ExponentialDisk   },
    { "miyamoto-nagai", MiyamotoNagaiDisk },
    END_MW_ENUM_ASSOCIATION
};

static int createDisk(lua_State* luaSt)
{
    Disk d = EMPTY_DISK;

    warn("Creating disk\n");

    d.type = checkEnum(luaSt, diskOptions, -1);
    pushDisk(luaSt, &d);
    return 1;
}

int getDiskT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, diskOptions, *(int*) v);
}

static int toStringDisk(lua_State* luaSt)
{
    Disk* d;
    char* str;

    d = checkDisk(luaSt, 1);
    str = showDisk(d);
    lua_pushstring(luaSt, str);
    free(str);

    return 1;
}

static const luaL_reg metaMethodsDisk[] =
{
    { "__tostring", toStringDisk },
    { NULL, NULL }
};

static const luaL_reg methodsDisk[] =
{
    { "create",   createDisk },
    { NULL, NULL }
};

static const Xet_reg_pre gettersDisk[] =
{
    { "type",         getDiskT,  offsetof(Disk, type) },
    { "mass" ,        getNumber, offsetof(Disk, mass) },
    { "scale_length", getNumber, offsetof(Disk, scale_length) },
    { "scale_height", getNumber, offsetof(Disk, scale_height) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersDisk[] =
{
    { "mass",         setNumber, offsetof(Disk, mass) },
    { "scale_length", setNumber, offsetof(Disk, scale_length) },
    { "scale_height", setNumber, offsetof(Disk, scale_height) },
    { NULL, NULL, 0 }
};

int registerDisk(lua_State* luaSt)
{
    return registerStruct(luaSt,
                          DISK_TYPE,
                          gettersDisk,
                          settersDisk,
                          metaMethodsDisk,
                          methodsDisk);
}

