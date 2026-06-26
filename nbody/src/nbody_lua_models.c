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
#include "nbody_dwarf_potential.h"

static const real pi = 3.1415926535;

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
            { "smalla", LUA_TNUMBER, NULL, TRUE,  &smalla, 1 },
            { "biga",   LUA_TNUMBER, NULL, TRUE,  &biga,   1 },
            { "Md",     LUA_TNUMBER, NULL, TRUE,  &Md,     1 },
            { "step",   LUA_TNUMBER, NULL, FALSE, &step,   1 },
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

real* nbCalculateEps2_NEW(const Dwarf* light_comp, const Dwarf* dark_comp, unsigned int lm_nbody, unsigned int nbody) //new softening length is dwarf specific and uses velocity dispersion approximation
{

    int dm_nbody = nbody - lm_nbody;
    if (dm_nbody == 0 || lm_nbody ==0)
        {
        Dwarf* comp = NULL;
        int bodies = 0;
        if (dm_nbody == 0) 
        {
            comp = light_comp;
            bodies = lm_nbody;
        }
        if (lm_nbody == 0) 
        {
            comp = dark_comp;
            bodies = dm_nbody;

        }
        // Separate code for 1-component systems
        real m = comp->mass / bodies;
        real rho_0 = get_density(comp, comp->scaleLength/10); //central density
        real d = 2* mw_pow(3*m/(4*pi*rho_0), 1.0/3.0);

        real radius = get_vel_disp_radius(comp);
        real m_encl = sqr(2*radius)*first_derivative(get_potential, 2*radius, comp)*-1;
        real v2 = m_encl / (2*radius);
        real r_strong = 2 * m / v2;

        real eps2 = r_strong * d;
        real eps = mw_sqrt(eps2);

        real* eps2_array = (real*)mwMalloc(sizeof(real) * 3);
        eps2_array[0] = eps2;
        eps2_array[1] = eps2;
        eps2_array[2] = eps2;
        mw_printf("Optimal Baryon Softening Length = %.15f kpc, Upper bound = %.15f kpc, Lower bound = %.15f kpc\n", eps, d, r_strong);
        return eps2_array;
        }


    // Average distance between stars
    
    real m_l = light_comp->mass / lm_nbody;
    real rho_0_l = get_density(light_comp, light_comp->scaleLength/10); //central density
    real d_l = 2* mw_pow(3*m_l/(4*pi*rho_0_l), 1.0/3.0);

    real m_d = dark_comp->mass / dm_nbody;
    real rho_0_d = get_density(dark_comp, dark_comp->scaleLength/10); //central density
    real d_d = 2* mw_pow(3*m_d/(4*pi*rho_0_d), 1.0/3.0);

    // Strong interaction radius
    real radius_l = get_vel_disp_radius(light_comp);
    real radius_d = get_vel_disp_radius(dark_comp);

    real m_encl_l = sqr(2*radius_l)*first_derivative(get_potential, 2*radius_l, light_comp)*-1 + sqr(radius_l)*first_derivative(get_potential, radius_l, dark_comp)*-1;
    real m_encl_d = sqr(2*radius_d)*first_derivative(get_potential, 2*radius_d, dark_comp)*-1 + sqr(radius_d)*first_derivative(get_potential, radius_d, light_comp)*-1;
    real v2_l = m_encl_l / (2*radius_l);
    real v2_d = m_encl_d / (2*radius_d);
    real r_strong_l = 2 * m_l / v2_l;
    real r_strong_d = 2 * m_d / v2_d;

    // Softening length
    // Suggestion for multidwarf: pass up both strong interaction radii and distance calculations and make the softening lengths in Lua (so you can do NxN)
    real cross_r_strong = (r_strong_l > r_strong_d) ? r_strong_l : r_strong_d;
    // Picks the larger strong interaction radius between the two (where the smaller particle will begin to strongly interact with the larger)
    real d_cross = (d_l < d_d) ? d_l : d_d;
    // Picks the smaller average distance (corresponding to the larger density)
    real eps2_l = r_strong_l * d_l;
    real eps_l = mw_sqrt(eps2_l);
    real eps2_d = r_strong_d * d_d;
    real eps_d = mw_sqrt(eps2_d);
    real eps2_cross = cross_r_strong * d_cross;
    real eps_cross = mw_sqrt(eps2_cross);

    mw_printf("Optimal Baryon Softening Length = %.15f kpc, Upper bound = %.15f kpc, Lower bound = %.15f kpc\n", eps_l, d_l, r_strong_l);
    mw_printf("Optimal Dark Matter Softening Length = %.15f kpc, Upper bound = %.15f kpc, Lower bound = %.15f kpc\n", eps_d, d_d, r_strong_d);
    mw_printf("Optimal Dark Matter-Baryon Softening Length = %.15f kpc, Upper bound = %.15f kpc, Lower bound = %.15f kpc\n", eps_cross, d_cross, cross_r_strong);
    real* eps2_array = (real*)mwMalloc(sizeof(real) * 3);
    eps2_array[0] = eps2_l;
    eps2_array[1] = eps2_cross;
    eps2_array[2] = eps2_d;
    return eps2_array;
}

static int luaCalculateEps2Dwarf(lua_State* luaSt) //read in params from lua to calc new softening length
{
    const Dwarf* model_1 = (const Dwarf*) mw_checknamedudata(luaSt, 1, DWARF_TYPE);
    const Dwarf* model_2 = (const Dwarf*) mw_checknamedudata(luaSt, 2, DWARF_TYPE);
    unsigned int lm_nbody = (unsigned int) luaL_checkinteger(luaSt, 3);
    unsigned int nbody = (unsigned int) luaL_checkinteger(luaSt, 4);
    
    real* array_ptr = nbCalculateEps2_NEW(model_1, model_2, lm_nbody, nbody);
    real eps2_l = array_ptr[0];
    real eps2_cross = array_ptr[1];
    real eps2_d = array_ptr[2];
    free(array_ptr);
    lua_pushnumber(luaSt, eps2_l);
    lua_pushnumber(luaSt, eps2_cross);
    lua_pushnumber(luaSt, eps2_d); // Return softening lengths for baryon-baryon, baryon-dark, and dark-dark interactions, respectively
    return 3;
}

__attribute__((unused)) static real nbCalculateEps2(real nbody, real a_b, real a_d, real M_b, real M_d) /* Eric's softening lenth, had some issues and is no longer used*/
{
    real beta = 1.0;                                  /** Tunable parameter for softening length **/
    real r_v = nbCalculateVirial(a_b, a_d, M_b, M_d); /** Calculate virial radius using formula for Henon length unit **/
    real eps = r_v * 0.98 * mw_pow(nbody, -0.26);     /** Optimal softening length pulled from Athanassoula et al. 1998 **/
    real eps2 = sqr(eps)/beta;
    if (eps2 <= REAL_EPSILON) {
        eps2 = REAL_EPSILON;
    }
    mw_printf("Optimal Softening Length = %.15f kpc\n", eps);
    return eps2;
}

static real nbCalculateEps2_OLD(real nbody, real a_b, real a_d, real M_b, real M_d) /** Old softening length formula from v1.76 and earlier **/
{
    real a = (M_b*a_b + M_d*a_d)/(M_b+M_d);
    real eps = a / 10.0 / mw_sqrt(nbody);
    real eps2 = sqr(eps);
    if (eps2 <= REAL_EPSILON) {
        eps2 = REAL_EPSILON;
    }
    mw_printf("Optimal Softening Length = %.15f kpc\n", eps);
    return eps2;
}

static int luaCalculateEps2_OLD(lua_State* luaSt) //read in params from lua to calc old softening length
{
    int nbody, arg_num;
    real a_b, a_d, M_b, M_d;

    arg_num = lua_gettop(luaSt);

    if (arg_num == 5)
    {
        nbody = (int) luaL_checkinteger(luaSt, 1);
        a_b = luaL_checknumber(luaSt, 2);
        a_d = luaL_checknumber(luaSt, 3);
        M_b = luaL_checknumber(luaSt, 4);
        M_d = luaL_checknumber(luaSt, 5);
        
    }
    else if (arg_num == 2) /** Single component only requires scale radius **/
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
    real eps2 = nbCalculateEps2_OLD(nbody, a_b, a_d, M_b, M_d);
    lua_pushnumber(luaSt, eps2);
    lua_pushnumber(luaSt, eps2);
    lua_pushnumber(luaSt, eps2); // Return same softening length for all interactions
    return 3;
}


static int luaReverseOrbit(lua_State* luaSt)
{
    mwvector finalPos = ZERO_VECTOR, finalVel = ZERO_VECTOR;
    static real dt = 0.0;
    static real tstop = 0.0;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;

    static const MWNamedArg argTable[] =
        {
            { "potential",  LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot,           1 },
            { "position",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos,           1 },
            { "velocity",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel,           1 },
            { "tstop",      LUA_TNUMBER,   NULL,           TRUE, &tstop,         1 },
            { "dt",         LUA_TNUMBER,   NULL,           TRUE, &dt,            1 },
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

static int luaReverseOrbit_LMC(lua_State* luaSt)
{
    mwvector finalPos = ZERO_VECTOR, finalVel = ZERO_VECTOR, LMCfinalPos = ZERO_VECTOR, LMCfinalVel = ZERO_VECTOR;
    static real dt = 0.0;
    static real tstop = 0.0;
    static real ftime = 0.0;
    static real LMCfunction = 1;
    static real LMCmass = 0.0;
    static real LMCscale = 0.0;
    static real LMCscale2 = 0.0;
    static real coulomb_log = 0.0;
    static mwbool LMCDynaFric = FALSE;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;
    static const mwvector* LMCpos = NULL;
    static const mwvector* LMCvel = NULL;

    static const MWNamedArg argTable[] =
        {
            { "potential",   LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot,         1 },
            { "position",    LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos,         1 },
            { "velocity",    LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel,         1 },
            { "LMCposition", LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &LMCpos,      1 },
            { "LMCvelocity", LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &LMCvel,      1 },
	        { "LMCfunction", LUA_TNUMBER,   NULL,           TRUE, &LMCfunction, 1 },
            { "LMCmass",     LUA_TNUMBER,   NULL,           TRUE, &LMCmass,     1 },
            { "LMCscale",    LUA_TNUMBER,   NULL,           TRUE, &LMCscale,    1 },
	        { "LMCscale2",   LUA_TNUMBER,   NULL,           TRUE, &LMCscale2,   1 },
            { "coulomb_log", LUA_TNUMBER,   NULL,           TRUE, &coulomb_log, 1 },
            { "LMCDynaFric", LUA_TBOOLEAN,  NULL,           TRUE, &LMCDynaFric, 1 },
            { "tstop",       LUA_TNUMBER,   NULL,           TRUE, &tstop,       1 },
            { "ftime",       LUA_TNUMBER,   NULL,           TRUE, &ftime,       1 },
            { "dt",          LUA_TNUMBER,   NULL,           TRUE, &dt,          1 },
            END_MW_NAMED_ARG
        };

    switch (lua_gettop(luaSt))
    {
        case 1:
            handleNamedArgumentTable(luaSt, argTable, 1);
            break;

        case 14:
            pot = checkPotential(luaSt, 1);
            pos = checkVector(luaSt, 2);
            vel = checkVector(luaSt, 3);
            LMCpos = checkVector(luaSt, 4);
            LMCvel = checkVector(luaSt, 5);
	        LMCfunction = luaL_checknumber(luaSt, 6);
            LMCmass = luaL_checknumber(luaSt, 7);
            LMCscale = luaL_checknumber(luaSt, 8);
	        LMCscale2 = luaL_checknumber(luaSt, 9);
            coulomb_log = luaL_checknumber(luaSt, 10);
            LMCDynaFric = luaL_checknumber(luaSt, 11);
            tstop = luaL_checknumber(luaSt, 12);
            ftime = luaL_checknumber(luaSt, 13);
            dt = luaL_checknumber(luaSt, 14);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 14 arguments");
    }

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    int lmcfunction = round(LMCfunction);
    nbReverseOrbit_LMC(&finalPos, &finalVel, &LMCfinalPos, &LMCfinalVel, pot, *pos, *vel, *LMCpos, *LMCvel, LMCDynaFric, ftime, tstop, dt, lmcfunction, LMCmass, LMCscale, LMCscale2, coulomb_log);
    pushVector(luaSt, finalPos);
    pushVector(luaSt, finalVel);
    pushVector(luaSt, LMCfinalPos);
    pushVector(luaSt, LMCfinalVel);

    return 4;
}

static int luaPrintReverseOrbit(lua_State* luaSt)
{
    mwvector finalPos = ZERO_VECTOR, finalVel = ZERO_VECTOR;
    static real dt = 0.0;
    static real tstop = 0.0;
    static real tstopf = 0.0;
    static Potential* pot = NULL;
    static const mwvector* pos = NULL;
    static const mwvector* vel = NULL;
    //static mwbool SecondDisk = FALSE;

    static const MWNamedArg argTable[] =
        {
            { "potential",  LUA_TUSERDATA, POTENTIAL_TYPE, TRUE, &pot,           1 },
            { "position",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &pos,           1 },
            { "velocity",   LUA_TUSERDATA, MWVECTOR_TYPE,  TRUE, &vel,           1 },
            { "tstop",      LUA_TNUMBER,   NULL,           TRUE, &tstop,         1 },
            { "tstopf",     LUA_TNUMBER,   NULL,           TRUE, &tstopf,        1 },
            { "dt",         LUA_TNUMBER,   NULL,           TRUE, &dt,            1 },
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
    lua_register(luaSt, "reverseOrbit_LMC", luaReverseOrbit_LMC);
    lua_register(luaSt, "PrintReverseOrbit", luaPrintReverseOrbit);
    lua_register(luaSt, "calculateEps2", luaCalculateEps2_OLD);
    lua_register(luaSt, "calculateEps2Dwarf", luaCalculateEps2Dwarf);
    lua_register(luaSt, "calculateTimestep", luaCalculateTimestep);
}
