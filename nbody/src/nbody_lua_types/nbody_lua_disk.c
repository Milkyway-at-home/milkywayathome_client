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

#include "nbody_types.h"
#include "nbody_show.h"
#include "nbody_lua_disk.h"
#include "milkyway_lua.h"
#include "milkyway_util.h"

Disk* checkDisk(lua_State* luaSt, int idx)
{
    return (Disk*) mw_checknamedudata(luaSt, idx, DISK_TYPE);
}

int pushDisk(lua_State* luaSt, const Disk* p)
{
    return pushType(luaSt, DISK_TYPE, sizeof(Disk), (void*) p);
}

static const MWEnumAssociation diskOptions[] =
{
    { "freeman",    FreemanDisk   },
    { "miyamoto-nagai", MiyamotoNagaiDisk },
    { "double-exponential", DoubleExponentialDisk },
    { "orbiting-bar", OrbitingBar}, //this is a time-dependent test potential
    { "none", NoDisk },
    END_MW_ENUM_ASSOCIATION
};

static int createDisk(lua_State* luaSt, const MWNamedArg* argTable, Disk* d)
{
    oneTableArgument(luaSt, argTable);
    pushDisk(luaSt, d);
    return 1;
}

static int createMiyamotoNagaiDisk(lua_State* luaSt)
{
    static Disk d = { MiyamotoNagaiDisk, 0.0, 0.0, 0.0 };

    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &d.mass        },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &d.scaleLength },
            { "scaleHeight", LUA_TNUMBER, NULL, TRUE, &d.scaleHeight },
            END_MW_NAMED_ARG
        };

    return createDisk(luaSt, argTable, &d);
}

static int createDoubleExponentialDisk(lua_State* luaSt)
{
    static Disk d = { DoubleExponentialDisk, 0.0, 0.0, 0.0 };

    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &d.mass        },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &d.scaleLength },
            { "scaleHeight", LUA_TNUMBER, NULL, TRUE, &d.scaleHeight },
            END_MW_NAMED_ARG
        };

    return createDisk(luaSt, argTable, &d);
}

static int createSech2ExponentialDisk(lua_State* luaSt)
{
    static Disk d = { Sech2ExponentialDisk, 0.0, 0.0, 0.0 };

    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &d.mass        },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &d.scaleLength },
            { "scaleHeight", LUA_TNUMBER, NULL, TRUE, &d.scaleHeight },
            END_MW_NAMED_ARG
        };

    return createDisk(luaSt, argTable, &d);
}

static int createFreemanDisk(lua_State* luaSt)
{
    static Disk d = { FreemanDisk, 0.0, 0.0, 0.0 };

    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &d.mass        },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &d.scaleLength },
            END_MW_NAMED_ARG
        };

    return createDisk(luaSt, argTable, &d);
}

static int createBar(lua_State* luaSt)
{
    static Disk d = { OrbitingBar, 0.0, 0.0, 0.0, 0.0 };

    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &d.mass        },
            { "scaleLength", LUA_TNUMBER, NULL, TRUE, &d.scaleLength },
            { "patternSpeed", LUA_TNUMBER, NULL, TRUE, &d.patternSpeed },
            { "startAngle", LUA_TNUMBER, NULL, TRUE, &d.startAngle },
            END_MW_NAMED_ARG
        };

    return createDisk(luaSt, argTable, &d);
}

static int createNoDisk(lua_State* luaSt)
{
    static Disk d = { NoDisk, 0.0, 0.0, 0.0 };

    static const MWNamedArg argTable[] =
        {
            { "mass",        LUA_TNUMBER, NULL, TRUE, &d.mass        },
            END_MW_NAMED_ARG
        };

    return createDisk(luaSt, argTable, &d);
}

int getDiskT(lua_State* luaSt, void* v)
{
    return pushEnum(luaSt, diskOptions, *(int*) v);
}

static int toStringDisk(lua_State* luaSt)
{
    return toStringType(luaSt, (StructShowFunc) showDisk, (LuaTypeCheckFunc) checkDisk);
}

static int eqDisk(lua_State* luaSt)
{
    lua_pushboolean(luaSt, equalDisk(checkDisk(luaSt, 1), checkDisk(luaSt, 2)));
    return 1;
}

int getDisk(lua_State* luaSt, void* v)
{
    pushDisk(luaSt, (Disk*) v);
    return 1;
}

int setDisk(lua_State* luaSt, void* v)
{
    *(Disk*) v = *checkDisk(luaSt, 3);
    return 0;
}

static const luaL_reg metaMethodsDisk[] =
{
    { "__tostring", toStringDisk },
    { "__eq",       eqDisk       },
    { NULL, NULL }
};

static const luaL_reg methodsDisk[] =
{
    { "miyamotoNagai", createMiyamotoNagaiDisk },
    { "freeman",   createFreemanDisk   },
    { "doubleExponential",   createDoubleExponentialDisk   },
    { "sech2Exponential",   createSech2ExponentialDisk   },
    { "orbitingBar", createBar  },
    { "none",   createNoDisk   },
    { NULL, NULL }
};

static const Xet_reg_pre gettersDisk[] =
{
    { "type",        getDiskT,  offsetof(Disk, type)        },
    { "mass" ,       getNumber, offsetof(Disk, mass)        },
    { "scaleLength", getNumber, offsetof(Disk, scaleLength) },
    { "scaleHeight", getNumber, offsetof(Disk, scaleHeight) },
    { "patternSpeed", getNumber, offsetof(Disk, patternSpeed) },
    { "startAngle", getNumber, offsetof(Disk, startAngle) },
    { NULL, NULL, 0 }
};

static const Xet_reg_pre settersDisk[] =
{
    { "mass",        setNumber, offsetof(Disk, mass)        },
    { "scaleLength", setNumber, offsetof(Disk, scaleLength) },
    { "scaleHeight", setNumber, offsetof(Disk, scaleHeight) },
    { "patternSpeed", getNumber, offsetof(Disk, patternSpeed) },
    { "startAngle", getNumber, offsetof(Disk, startAngle) },
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

int registerDiskKinds(lua_State* luaSt)
{
    int table;

    lua_newtable(luaSt);
    table = lua_gettop(luaSt);

    setModelTableItem(luaSt, table, createMiyamotoNagaiDisk, "miyamotoNagai");
    setModelTableItem(luaSt, table, createFreemanDisk, "freeman");
    setModelTableItem(luaSt, table, createDoubleExponentialDisk, "doubleExponential");
    setModelTableItem(luaSt, table, createSech2ExponentialDisk, "sech2Exponential");
    setModelTableItem(luaSt, table, createBar, "orbitingBar");
    setModelTableItem(luaSt, table, createNoDisk, "none");

    lua_setglobal(luaSt, "diskModels");

    return 0;
}

