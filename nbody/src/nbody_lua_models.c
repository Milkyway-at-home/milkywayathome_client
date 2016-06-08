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
#include "nbody_nfw.h"
#include "nbody_hernq.h"
#include "nbody_isotropic.h"
#include "nbody_lua_models.h"
#include "nbody_check_params.h"
#include "nbody_defaults.h"


/* For using a combination of light and dark models to generate timestep */
static real plummerTimestepIntegral(real smalla, real biga, real Md, real step)
{
    /* Calculate the enclosed mass of the big sphere within the little sphere's scale length */
    real encMass, val, r;

    encMass = 0.0;
    for (r = 0.0; r <= smalla; r += step)
    {
        val = sqr(r) / mw_pow(sqr(r) + sqr(biga), 2.5);
        encMass += val * step;
    }
    encMass *= 3.0 * Md * sqr(biga);

    return encMass;
}

static int luaPlummerTimestepIntegral(lua_State* luaSt)
{
    int nArgs;
    static real smalla = 0.0, biga = 0.0, Md = 0.0, encMass = 0.0;
    static real step = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "smalla", LUA_TNUMBER, NULL, TRUE,  &smalla },
            { "biga",   LUA_TNUMBER, NULL, TRUE,  &biga   },
            { "Md",     LUA_TNUMBER, NULL, TRUE,  &Md     },
            { "step",   LUA_TNUMBER, NULL, FALSE, &step   },
            END_MW_NAMED_ARG
        };

    step = 1.0e-5;  /* Reset default step */
    nArgs = lua_gettop(luaSt);
    if (nArgs == 1 && lua_istable(luaSt, 1))
    {
        handleNamedArgumentTable(luaSt, argTable, 1);
    }
    else if (nArgs == 3 || nArgs == 4)
    {
        smalla = luaL_checknumber(luaSt, 1);
        biga = luaL_checknumber(luaSt, 2);
        Md = luaL_checknumber(luaSt, 3);
        step = luaL_optnumber(luaSt, 4, step);
    }
    else
    {
        return luaL_argerror(luaSt, 1, "Expected 1, 3 or 4 arguments");
    }

    /* Make sure the bounds / step are OK so that this integral will be sure to complete */
    if (mwCheckNormalPosNum(smalla))
        return luaL_argerror(luaSt, 1, "Invalid small radius");
    if (mwCheckNormalPosNum(biga))
        return luaL_argerror(luaSt, 2, "Invalid big radius");
    if (mwCheckNormalPosNumEps(step))
        return luaL_argerror(luaSt, 4, "Invalid step argument");

    encMass = plummerTimestepIntegral(smalla, biga, Md, step);
    lua_pushnumber(luaSt, encMass);

    return 1;
}

void registerPredefinedModelGenerators(lua_State* luaSt)
{
    int table;

    registerGeneratePlummer(luaSt);
    registerGenerateNFW(luaSt);
    registerGenerateHernq(luaSt);
    registerGenerateIsotropic(luaSt);

    /* Create a table of predefined models, so we can use them like
     * predefinedModels.plummer() etc. */
    lua_newtable(luaSt);
    table = lua_gettop(luaSt);

    setModelTableItem(luaSt, table, nbGeneratePlummer, "plummer");
    setModelTableItem(luaSt, table, nbGenerateNFW, "nfw");
    setModelTableItem(luaSt, table, nbGenerateHernq, "hernq");
    setModelTableItem(luaSt, table, nbGenerateIsotropic, "isotropic");
    /*
      setModelTableItem(luaSt, table, generateKing, "king");
      setModelTableItem(luaSt, table, generateDehnen, "dehnen");
    */

    lua_setglobal(luaSt, "predefinedModels");
}

static real nbCalculateTimestep(real mass, real r0)
{
    return sqr(1.0/10.0) * mw_sqrt((PI_4_3 * cube(r0)) / mass);
}

static int luaCalculateTimestep(lua_State* luaSt)
{
    real mass, r0;

    if (lua_gettop(luaSt) != 2)
        return luaL_argerror(luaSt, 0, "Expected 2 arguments");

    mass = luaL_checknumber(luaSt, 1);
    r0 = luaL_checknumber(luaSt, 2);

    lua_pushnumber(luaSt, nbCalculateTimestep(mass, r0));
    return 1;
}

static real nbCalculateEps2(real nbody, real r0)
{
    real eps = r0 / (10.0 * mw_sqrt(nbody));
    return sqr(eps);
}

static int luaCalculateEps2(lua_State* luaSt)
{
    int nbody;
    real r0;

    if (lua_gettop(luaSt) != 2)
        return luaL_argerror(luaSt, 0, "Expected 2 arguments");

    nbody = (int) luaL_checkinteger(luaSt, 1);
    r0 = luaL_checknumber(luaSt, 2);

    lua_pushnumber(luaSt, nbCalculateEps2((real) nbody, r0));

    return 1;
}

static int luaReverseOrbit(lua_State* luaSt)
{
    mwvector finalPos, finalVel;
    static real dt = 0.0;
    static real tstop = 0.0;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;

    static const MWNamedArg argTable[] =
        {
            { "potential", LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot   },
            { "position",  LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos   },
            { "velocity",  LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel   },
            { "tstop",     LUA_TNUMBER,   NULL,           TRUE, &tstop },
            { "dt",        LUA_TNUMBER,   NULL,           TRUE, &dt    },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 5:
            pot = checkPotential(luaSt, 1);
            pos = checkVector(luaSt, 2);
            vel = checkVector(luaSt, 3);
            tstop = luaL_checknumber(luaSt, 4);
            dt = luaL_checknumber(luaSt, 5);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 5 arguments");
    }

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    nbReverseOrbit(&finalPos, &finalVel, pot, *pos, *vel, tstop, dt);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);

    return 2;
}


static int luaPrintReverseOrbit(lua_State* luaSt)
{
    mwvector finalPos, finalVel;
    static real dt = 0.0;
    static real tstop = 0.0;
    static real tstopf = 0.0;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;

    static const MWNamedArg argTable[] =
        {
            { "potential", LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot    },
            { "position",  LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos    },
            { "velocity",  LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel    },
            { "tstop",     LUA_TNUMBER,   NULL,           TRUE, &tstop  },
            { "tstopf",    LUA_TNUMBER,   NULL,           TRUE, &tstopf },
            { "dt",        LUA_TNUMBER,   NULL,           TRUE, &dt     },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 6:
            pot = checkPotential(luaSt, 1);
            pos = checkVector(luaSt, 2);
            vel = checkVector(luaSt, 3);
            tstop = luaL_checknumber(luaSt, 4);
            tstopf = luaL_checknumber(luaSt, 5);
            dt = luaL_checknumber(luaSt, 6);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 6 arguments");
    }

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    nbPrintReverseOrbit(&finalPos, &finalVel, pot, *pos, *vel, tstop, tstopf, dt);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);

    return 2;
}

void registerModelFunctions(lua_State* luaSt)
{
    lua_register(luaSt, "plummerTimestepIntegral", luaPlummerTimestepIntegral);
    lua_register(luaSt, "reverseOrbit", luaReverseOrbit);
    lua_register(luaSt, "PrintReverseOrbit", luaPrintReverseOrbit);
    lua_register(luaSt, "calculateEps2", luaCalculateEps2);
    lua_register(luaSt, "calculateTimestep", luaCalculateTimestep);
}

