/*
 *  Copyright (c) 2011 Matthew Arsenault
 *  Copyright (c) 2016-2018 Siddhartha Shelton
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
#include "nbody_virial.h"
#include "nbody_plummer.h"
#include "nbody_nfw.h"
#include "nbody_hernq.h"
#include "nbody_isotropic.h"
#include "nbody_mixeddwarf.h"
#include "nbody_manual_bodies.h"
#include "nbody_lua_models.h"
#include "nbody_check_params.h"
#include "nbody_defaults.h"
#include "nbody_potential_types.h"
#include "nbody_lua_dwarf.h"
#include "nbody_autodiff.h"

/* For using a combination of light and dark models to generate timestep */
static real_0 plummerTimestepIntegral(real_0 smalla, real_0 biga, real_0 Md, real_0 step)
{
    /* Calculate the enclosed mass of the big sphere within the little sphere's scale length */
    real_0 encMass, val, r;

    encMass = 0.0;
    for (r = 0.0; r <= smalla; r += step)
    {
        val = sqr_0(r) / mw_pow_0(sqr_0(r) + sqr_0(biga), 2.5);
        encMass += val * step;
    }
    encMass *= 3.0 * Md * sqr_0(biga);

    return encMass;
}

static int luaPlummerTimestepIntegral(lua_State* luaSt)
{
    int nArgs;
    static real_0 smalla = 0.0, biga = 0.0, Md = 0.0, encMass = 0.0;
    static real_0 step = 0.0;

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
    registerGenerateMixedDwarf(luaSt);
    registerGenerateManualBodies(luaSt);
    
    /* Create a table of predefined models, so we can use them like
     * predefinedModels.plummer() etc. */
    lua_newtable(luaSt);
    table = lua_gettop(luaSt);

    setModelTableItem(luaSt, table, nbGeneratePlummer, "plummer");
    setModelTableItem(luaSt, table, nbGenerateNFW, "nfw");
    setModelTableItem(luaSt, table, nbGenerateHernq, "hernq");
    setModelTableItem(luaSt, table, nbGenerateIsotropic, "isotropic");
    setModelTableItem(luaSt, table, nbGenerateMixedDwarf, "mixeddwarf");
    setModelTableItem(luaSt, table, nbGenerateManualBodies, "manual_bodies");
    
    /*
      setModelTableItem(luaSt, table, generateKing, "king");
      setModelTableItem(luaSt, table, generateDehnen, "dehnen");
    */

    lua_setglobal(luaSt, "predefinedModels");
}

static real_0 nbCalculateTimestep(real_0 mass, real_0 r0)
{
    return sqr_0(1.0/10.0) * mw_sqrt_0((PI_4_3 * cube_0(r0)) / mass);
}

static int luaCalculateTimestep(lua_State* luaSt)
{
    real_0 mass, r0;

    if (lua_gettop(luaSt) != 2)
        return luaL_argerror(luaSt, 0, "Expected 2 arguments");

    mass = luaL_checknumber(luaSt, 1);
    r0 = luaL_checknumber(luaSt, 2);

    lua_pushnumber(luaSt, nbCalculateTimestep(mass, r0));
    return 1;
}

static real_0 nbCalculateEps2(real_0 nbody, real_0 a_b, real_0 a_d, real_0 M_b, real_0 M_d)
{
    real_0 beta = 1.0;                                  /** Tunable parameter for softening length **/
    real_0 r_v = nbCalculateVirial(a_b, a_d, M_b, M_d); /** Calculate virial radius using formula for Henon length unit **/
    real_0 eps = r_v * 0.98 * mw_pow_0(nbody, -0.26);     /** Optimal softening length pulled from Athanassoula et al. 1998 **/
    real_0 eps2 = sqr_0(eps)/beta;
    if (eps2 <= REAL_EPSILON) {
        eps2 = REAL_EPSILON;
    }
    //mw_printf("Optimal Softening Length = %.15f kpc\n", eps);
    return eps2;
}

static int luaCalculateEps2(lua_State* luaSt)
{
    int nbody, arg_num;
    real_0 r0, a_b, a_d, M_b, M_d;

    arg_num = lua_gettop(luaSt);

    if (arg_num == 5)
    {
        nbody = (int) luaL_checkinteger(luaSt, 1);
        a_b = luaL_checknumber(luaSt, 2);
        a_d = luaL_checknumber(luaSt, 3);
        M_b = luaL_checknumber(luaSt, 4);
        M_d = luaL_checknumber(luaSt, 5);
    }
    else if (arg_num == 2) /** Single component only requires scale radius. **/
    {
        nbody = (int) luaL_checkinteger(luaSt, 1);
        a_b = luaL_checknumber(luaSt, 2);
        a_d = 1.0; /** can be anything but zero **/
        M_b = 1.0; /** can be anything but zero **/
        M_d = 0.0; /** must be zero **/
    }
    else
    {
        return luaL_argerror(luaSt, 0, "Expected 2 or 5 arguments");
    }

    lua_pushnumber(luaSt, nbCalculateEps2((real_0) nbody, a_b, a_d, M_b, M_d));

    return 1;
}

static int luaReverseOrbit(lua_State* luaSt)
{
    mwvector finalPos, finalVel;
    static real_0 dt = 0.0;
    static real_0 tstop = 0.0;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;
    static real_0 sun_dist = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "potential",  LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot           },
            { "position",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos           },
            { "velocity",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel           },
            { "tstop",      LUA_TNUMBER,   NULL,           TRUE, &tstop         },
            { "dt",         LUA_TNUMBER,   NULL,           TRUE, &dt            },
            { "sunGCDist",  LUA_TNUMBER,   NULL,           TRUE, &sun_dist      },
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
            dt = luaL_checknumber(luaSt, 5);
            sun_dist = luaL_checknumber(luaSt, 6);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 6 arguments");
    }

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    //mw_printf("POS = [%.15f, %.15f, %.15f]\n", showRealValue(&pos->x), showRealValue(&pos->y), showRealValue(&pos->z));
    nbReverseOrbit(&finalPos, &finalVel, pot, pos, vel, tstop, dt, sun_dist);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);

    return 2;
}

static int luaReverseOrbit_LMC(lua_State* luaSt)
{
    mwvector finalPos, finalVel, LMCfinalPos, LMCfinalVel;
    static real_0 dt = 0.0;
    static real_0 tstop = 0.0;
    static real_0 ftime = 0.0;
    static real_0 LMCmass = 0.0;
    static real_0 LMCscale = 0.0;
    static real_0 coulomb_log = 0.0;
    static mwbool LMCDynaFric = FALSE;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;
    static const mwvector* LMCpos = NULL;
    static const mwvector* LMCvel = NULL;
    static real_0 sun_dist = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "potential",   LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot         },
            { "position",    LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos         },
            { "velocity",    LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel         },
            { "LMCposition", LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &LMCpos      },
            { "LMCvelocity", LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &LMCvel      },
            { "LMCmass",     LUA_TNUMBER,   NULL,           TRUE, &LMCmass     },
            { "LMCscale",    LUA_TNUMBER,   NULL,           TRUE, &LMCscale    },
            { "coulomb_log", LUA_TNUMBER,   NULL,           TRUE, &coulomb_log },
            { "LMCDynaFric", LUA_TBOOLEAN,  NULL,           TRUE, &LMCDynaFric },
            { "tstop",       LUA_TNUMBER,   NULL,           TRUE, &tstop       },
            { "ftime",       LUA_TNUMBER,   NULL,           TRUE, &ftime       },
            { "dt",          LUA_TNUMBER,   NULL,           TRUE, &dt          },
            { "sunGCDist",   LUA_TNUMBER,   NULL,           TRUE, &sun_dist    },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 13:
            pot = checkPotential(luaSt, 1);
            pos = checkVector(luaSt, 2);
            vel = checkVector(luaSt, 3);
            LMCpos = checkVector(luaSt, 4);
            LMCvel = checkVector(luaSt, 5);
            LMCmass = luaL_checknumber(luaSt, 6);
            LMCscale = luaL_checknumber(luaSt, 7);
            coulomb_log = luaL_checknumber(luaSt, 8);
            LMCDynaFric = luaL_checknumber(luaSt, 9);
            tstop = luaL_checknumber(luaSt, 10);
            ftime = luaL_checknumber(luaSt, 11);
            dt = luaL_checknumber(luaSt, 12);
            sun_dist = luaL_checknumber(luaSt, 13);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 13 arguments");
    }

    real LMCmass_var = mw_real_var(LMCmass, LMC_MASS_POS);
    real LMCscale_var = mw_real_var(LMCscale, LMC_RADIUS_POS);

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    nbReverseOrbit_LMC(&finalPos, &finalVel, &LMCfinalPos, &LMCfinalVel, pot, pos, vel, LMCpos, LMCvel, LMCDynaFric, ftime, tstop, dt, &LMCmass_var, &LMCscale_var, sun_dist, coulomb_log);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);
    pushVector(luaSt, LMCfinalPos);
    pushVector(luaSt, LMCfinalVel);

    return 4;
}

static int luaPrintReverseOrbit(lua_State* luaSt)
{
    mwvector finalPos, finalVel;
    static real_0 dt = 0.0;
    static real_0 tstop = 0.0;
    static real_0 tstopf = 0.0;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;
    //static mwbool SecondDisk = FALSE;

    static const MWNamedArg argTable[] =
        {
            { "potential",  LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot           },
            { "position",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos           },
            { "velocity",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel           },
            { "tstop",      LUA_TNUMBER,   NULL,           TRUE, &tstop         },
            { "tstopf",     LUA_TNUMBER,   NULL,           TRUE, &tstopf        },
            { "dt",         LUA_TNUMBER,   NULL,           TRUE, &dt            },
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

    nbPrintReverseOrbit(&finalPos, &finalVel, pot, pos, vel, tstop, tstopf, dt);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);

    return 2;
}

static int luaPrintReverseOrbit_LMC(lua_State* luaSt)
{
    mwvector finalPos, finalVel, LMCfinalPos, LMCfinalVel;
    static real_0 coulomb_log = 0.0;
    static real_0 dt = 0.0;
    static real_0 tstop = 0.0;
    static real_0 tstopf = 0.0;
    static real_0 LMCmass = 0.0;
    static real_0 LMCscale = 0.0;
    static mwbool LMCDynaFric = FALSE;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;
    static const mwvector* LMCpos = NULL;
    static const mwvector* LMCvel = NULL;

    static const MWNamedArg argTable[] =
        {
            { "potential",   LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot         },
            { "position",    LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos         },
            { "velocity",    LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel         },
            { "LMCposition", LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &LMCpos      },
            { "LMCvelocity", LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &LMCvel      },
            { "LMCmass",     LUA_TNUMBER,   NULL,           TRUE, &LMCmass     },
            { "LMCscale",    LUA_TNUMBER,   NULL,           TRUE, &LMCscale    },
            { "LMCDynaFric", LUA_TBOOLEAN,  NULL,           TRUE, &LMCDynaFric },
            { "tstop",       LUA_TNUMBER,   NULL,           TRUE, &tstop       },
            { "tstopf",      LUA_TNUMBER,   NULL,           TRUE, &tstopf      },
            { "dt",          LUA_TNUMBER,   NULL,           TRUE, &dt          },
            { "coulomb_log", LUA_TNUMBER,   NULL,           TRUE, &coulomb_log },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 12:
            pot = checkPotential(luaSt, 1);
            pos = checkVector(luaSt, 2);
            vel = checkVector(luaSt, 3);
            LMCpos = checkVector(luaSt, 4);
            LMCvel = checkVector(luaSt, 5);
            LMCmass = luaL_checknumber(luaSt, 6);
            LMCscale = luaL_checknumber(luaSt, 7);
            LMCDynaFric = luaL_checknumber(luaSt, 8);
            tstop = luaL_checknumber(luaSt, 9);
            tstopf = luaL_checknumber(luaSt, 10);
            dt = luaL_checknumber(luaSt, 11);
            coulomb_log = luaL_checknumber(luaSt, 12);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 12 arguments");
    }

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    real LMCmass_var = mw_real_var(LMCmass, LMC_MASS_POS);
    real LMCscale_var = mw_real_var(LMCscale, LMC_RADIUS_POS);

    nbPrintReverseOrbit_LMC(&finalPos, &finalVel, &LMCfinalPos, &LMCfinalVel, pot, pos, vel, LMCpos, LMCvel, LMCDynaFric, tstop, tstopf, dt, &LMCmass_var, &LMCscale_var, coulomb_log);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);
    pushVector(luaSt, LMCfinalPos);
    pushVector(luaSt, LMCfinalVel);

    return 4;
}

void registerModelFunctions(lua_State* luaSt)
{
    lua_register(luaSt, "plummerTimestepIntegral", luaPlummerTimestepIntegral);
    lua_register(luaSt, "reverseOrbit", luaReverseOrbit);
    lua_register(luaSt, "reverseOrbit_LMC", luaReverseOrbit_LMC);
    lua_register(luaSt, "PrintReverseOrbit", luaPrintReverseOrbit);
    lua_register(luaSt, "PrintReverseOrbit_LMC", luaPrintReverseOrbit_LMC);
    lua_register(luaSt, "calculateEps2", luaCalculateEps2);
    lua_register(luaSt, "calculateTimestep", luaCalculateTimestep);
}
