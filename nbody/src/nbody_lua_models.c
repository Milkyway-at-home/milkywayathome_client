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

static real nbCalculateEps2_NEW(const Dwarf* light_comp, unsigned int nbody) //new softening length is dwarf specific and needs velocity dispersion
{
    // Average distance between stars
    real m = light_comp->mass / nbody;
    real rho_0 = get_density(light_comp, 0.00001);
    real d = 2* mw_pow(3*m/(4*pi*rho_0), 1.0/3.0);

    // Strong interaction radius
    real v_disp = get_vel_disp(light_comp);
    real r_strong = 2 * m / v_disp; 

    // Softening length
    real logeps = 0.5 * (mw_log10(d) + mw_log10(r_strong));
    real eps = mw_pow(10, logeps);
    real eps2 = sqr(eps);
    mw_printf("Optimal Softening Length = %.15f kpc, Upper bound = %.15f kpc, Lower bound = %.15f kpc\n", eps, d, r_strong);
    return eps2;
}

static int luaCalculateEps2Dwarf(lua_State* luaSt) //read in params from lua to calc new softening length
{
    const Dwarf* model = (const Dwarf*) mw_checknamedudata(luaSt, 1, DWARF_TYPE);
    unsigned int nbody = (unsigned int) luaL_checkinteger(luaSt, 2);
    
    lua_pushnumber(luaSt, nbCalculateEps2_NEW(model, nbody));
    return 1;
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
    lua_pushnumber(luaSt, nbCalculateEps2_OLD((real) nbody, a_b, a_d, M_b, M_d));
    return 1;
}


static int luaReverseOrbitS(lua_State* luaSt)
{       
    size_t lenPos = lua_objlen(luaSt, 2); 
    size_t lenVel = lua_objlen(luaSt, 3);

    int arg_num;
    Potential* pot;
    mwvector pos[lenPos];
    mwvector vel[lenPos];
    mwvector finalPos[lenPos];
    mwvector finalVel[lenPos];
    real tstop, dt;
    real masses[lenPos];

    arg_num = lua_gettop(luaSt);

    if (arg_num != 6)
    {
        return luaL_argerror(luaSt, 0, "Expected 6 arguments");
    }


    if (lenPos != lenVel) {
        return luaL_error(luaSt, "Lengths of tables do not match");
    }

    pot = (Potential*) luaL_checkudata(luaSt, 1, POTENTIAL_TYPE);
    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    for (int i = 1; i <= lenPos; ++i) 
    {
        lua_rawgeti(luaSt, 2, i);  
        pos[i-1] = *(mwvector*) checkVector(luaSt, -1); 
        lua_pop(luaSt, 1);

        lua_rawgeti(luaSt, 3, i);  
        vel[i-1] = *(mwvector*) checkVector(luaSt, -1); 
        lua_pop(luaSt, 1);

        lua_rawgeti(luaSt, 6, i);
        masses[i-1] = luaL_checknumber(luaSt, -1);
        lua_pop(luaSt, 1);
    }

    tstop = luaL_checknumber(luaSt, 4);
    dt = luaL_checknumber(luaSt, 5);
    
    nbReverseOrbitS(finalPos, finalVel, pot, pos, vel, lenPos, tstop, dt, masses);

    pushVectorTable(luaSt, finalPos, lenPos);
    pushVectorTable(luaSt, finalVel, lenPos);

    return 2;
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

static int luaReverseOrbitS_LMC(lua_State* luaSt)
{
    size_t lenPos = lua_objlen(luaSt, 2); 
    size_t lenVel = lua_objlen(luaSt, 3);
    
    int arg_num;
    Potential* pot;
    mwvector pos[lenPos];
    mwvector vel[lenPos];
    mwvector finalPos[lenPos];
    mwvector finalVel[lenPos];
    mwvector* LMCpos; 
    mwvector* LMCvel;
    mwvector LMCfinalPos, LMCfinalVel;
    mwbool LMCDynaFric = FALSE;
    real LMCmass, LMCscale, coulomb_log, tstop, dt, ftime;
    real masses[lenPos];
    real rscales[lenPos];

    arg_num = lua_gettop(luaSt);

    if (arg_num != 14)
    {
        return luaL_argerror(luaSt, 0, "Expected 14 arguments");
    }


    if (lenPos != lenVel) {
        return luaL_error(luaSt, "Lengths of tables do not match");
    }

    for (int i = 1; i <= lenPos; ++i) 
    {
        lua_rawgeti(luaSt, 2, i);  
        pos[i-1] = *(mwvector*) checkVector(luaSt, -1); 
        lua_pop(luaSt, 1);

        lua_rawgeti(luaSt, 3, i);  
        vel[i-1] = *(mwvector*) checkVector(luaSt, -1); 
        lua_pop(luaSt, 1);

        lua_rawgeti(luaSt, 13, i);
        masses[i-1] = luaL_checknumber(luaSt, -1);
        lua_pop(luaSt, 1);

        lua_rawgeti(luaSt, 14, i);
        rscales[i-1] = luaL_checknumber(luaSt, -1);
        lua_pop(luaSt, 1);
        
    }

    pot = (Potential*) luaL_checkudata(luaSt, 1, POTENTIAL_TYPE);
    LMCpos = checkVector(luaSt, 4);
    LMCvel = checkVector(luaSt, 5);
    LMCmass = luaL_checknumber(luaSt, 6);
    LMCscale = luaL_checknumber(luaSt, 7);
    LMCDynaFric = luaL_checknumber(luaSt, 8);
    coulomb_log = luaL_checknumber(luaSt, 9);
    tstop = luaL_checknumber(luaSt, 10);
    ftime = luaL_checknumber(luaSt, 11);
    dt = luaL_checknumber(luaSt, 12);


    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    nbReverseOrbitS_LMC(finalPos, finalVel, &LMCfinalPos, &LMCfinalVel, pot, pos, vel, lenPos, *LMCpos, *LMCvel, LMCDynaFric, ftime, tstop, dt, LMCmass, LMCscale, coulomb_log, masses, rscales);
    pushVectorTable(luaSt, finalPos, lenPos);
    pushVectorTable(luaSt, finalVel, lenPos);
    pushVector(luaSt, LMCfinalPos);
    pushVector(luaSt, LMCfinalVel);

    return 4;
}

static int luaReverseOrbit_LMC(lua_State* luaSt)
{
    mwvector finalPos = ZERO_VECTOR, finalVel = ZERO_VECTOR, LMCfinalPos = ZERO_VECTOR, LMCfinalVel = ZERO_VECTOR;
    static real dt = 0.0;
    static real tstop = 0.0;
    static real ftime = 0.0;
    static real LMCmass = 0.0;
    static real LMCscale = 0.0;
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
            { "LMCmass",     LUA_TNUMBER,   NULL,           TRUE, &LMCmass,     1 },
            { "LMCscale",    LUA_TNUMBER,   NULL,           TRUE, &LMCscale,    1 },
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

        case 12:
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
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 1 or 12 arguments");
    }

    /* Make sure precalculated constants ready for use */
    if (checkPotentialConstants(pot))
        luaL_error(luaSt, "Error with potential");

    nbReverseOrbit_LMC(&finalPos, &finalVel, &LMCfinalPos, &LMCfinalVel, pot, *pos, *vel, *LMCpos, *LMCvel, LMCDynaFric, ftime, tstop, dt, LMCmass, LMCscale, coulomb_log);
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
    lua_register(luaSt, "reverseOrbitS", luaReverseOrbitS);
    lua_register(luaSt, "reverseOrbit", luaReverseOrbit);    
    lua_register(luaSt, "reverseOrbitS_LMC", luaReverseOrbitS_LMC);
    lua_register(luaSt, "reverseOrbit_LMC", luaReverseOrbit_LMC);
    lua_register(luaSt, "PrintReverseOrbit", luaPrintReverseOrbit);
    lua_register(luaSt, "calculateEps2", luaCalculateEps2_OLD);
    lua_register(luaSt, "calculateEps2Dwarf", luaCalculateEps2Dwarf);
    lua_register(luaSt, "calculateTimestep", luaCalculateTimestep);
}
