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
#include "milkyway_util.h"
#include "nbody_orbit_integrator.h"
#include "nbody_coordinates.h"

#include "nbody_lua.h"
#include "nbody_lua_types.h"
#include "milkyway_lua.h"
#include "nbody_plummer.h"
#include "nbody_lua_util.h"
#include "nbody_check_params.h"
#include "nbody_defaults.h"


static int luaLbrToCartesian(lua_State* luaSt)
{
    mwbool useRadians = FALSE, useGalacticCoordinates = FALSE;
    real sunGCDist = DEFAULT_SUN_GC_DISTANCE;
    const NBodyCtx* ctx = NULL;
    mwvector v;

    if (lua_gettop(luaSt) > 4)
        luaL_argerror(luaSt, 4, "Expected 1, 2, 3 or 4 arguments");

    ctx = checkNBodyCtx(luaSt, 1);
    v = *checkVector(luaSt, 2);

    /* ctx = toNBodyCtx(luaSt, 2);
       sunGCDist = ctx != NULL ? ctx->sunGCDist : luaL_optnumber(luaSt, 2, DEFAULT_SUN_GC_DISTANCE);
    */

    sunGCDist = ctx->sunGCDist;

    useGalacticCoordinates = mw_lua_optboolean(luaSt, 3, FALSE);
    useRadians = mw_lua_optboolean(luaSt, 4, FALSE);

    if (!useGalacticCoordinates)
        v = useRadians ? lbrToCartesian_rad(v, sunGCDist) : lbrToCartesian(v, sunGCDist);

    pushVector(luaSt, v);

    return 1;
}


static const real solarMassesPerMassUnit = 222288.47;

static int luaSolarMassToMassUnit(lua_State* luaSt)
{
    real mSolar;
    real m;

    if (lua_gettop(luaSt) != 1)
    {
        return luaL_argerror(luaSt, 0, "Expected 1 argument");
    }

    mSolar = luaL_checknumber(luaSt, 1);
    m = mSolar / solarMassesPerMassUnit;
    lua_pushnumber(luaSt, m);

    return 1;
}

static int luaMassUnitToSolarMass(lua_State* luaSt)
{
    real mSolar;
    real m;

    if (lua_gettop(luaSt) != 1)
    {
        return luaL_argerror(luaSt, 0, "Expected 1 argument");
    }

    m = luaL_checknumber(luaSt, 1);
    mSolar = solarMassesPerMassUnit * m;
    lua_pushnumber(luaSt, mSolar);

    return 1;
}

static int luaLightyearToKiloparsec(lua_State* luaSt)
{
    real kpc;
    real ly;

    if (lua_gettop(luaSt) != 1)
    {
        return luaL_argerror(luaSt, 0, "Expected 1 argument");
    }

    ly = luaL_checknumber(luaSt, 1);
    kpc = 0.000306599 * ly;
    lua_pushnumber(luaSt, kpc);

    return 1;

}

static int luaKiloparsecToLightyear(lua_State* luaSt)
{
    real kpc;
    real ly;

    if (lua_gettop(luaSt) != 1)
    {
        return luaL_argerror(luaSt, 0, "Expected 1 argument");
    }

    kpc = luaL_checknumber(luaSt, 1);
    ly = 3271.59 * kpc;
    lua_pushnumber(luaSt, ly);

    return 1;

}

void nbRegisterUtilityFunctions(lua_State* luaSt)
{
    lua_register(luaSt, "lbrToCartesian", luaLbrToCartesian);

    /* Unit conversions */
    lua_register(luaSt, "solarMassToMassUnit", luaSolarMassToMassUnit);
    lua_register(luaSt, "massUnitToSolarMass", luaMassUnitToSolarMass);
    lua_register(luaSt, "lightyearToKiloparsec", luaLightyearToKiloparsec);
    lua_register(luaSt, "kiloparsecToLightyear", luaKiloparsecToLightyear);
}


