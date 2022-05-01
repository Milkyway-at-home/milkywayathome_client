/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
   Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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

#include "nbody_priv.h"
#include "milkyway_util.h"
#include "milkyway_lua.h"
#include "nbody_lua_types.h"
#include "nbody_nfw.h"

//FIXME: THIS CODE DOES NOT PROPAGATE DERIVAIVE INFORMATION! MUST BE REWORKED BEFORE RUNNING WITH AUTODIFF!

static real_0 nfwMassInsideRadius(real_0 radius, real_0 rho_0, real_0 R_S)
{
    //Returns that mass inside a certain radius

    //Leading Constant
    real_0 mass = 4.0 * M_PI * rho_0 * cube_0(radius);
    //Integration terms
    mass *= ( mw_log_0((R_S + radius) / R_S) - radius / (R_S + radius));
    return mass;

}

static real_0 nfwNextRadius(real_0 start_radius, real_0 goal_mass, real_0 rho_0, real_0 R_S)
{
    // This is a scary function which returns the next radius limit
    // Do not use this code on client computers. It is slow inneficient
    // and well just dirty.
    real_0 test_radius;

    for (test_radius = start_radius;; test_radius += (real_0)0.0001)
    {
        if (nfwMassInsideRadius(test_radius, rho_0, R_S) >= goal_mass)
        {
            return test_radius;
        }
    }
}

/* nfwPickShell: pick a random point on a sphere of specified radius. */
static mwvector nfwPickShell(dsfmt_t* dsfmtState, real_0 rad)
{
    real_0 rsq;
    real_0 rsc;
    mwvector vec;
    real tmp;

    do                      /* pick point in NDIM-space */
    {
        vec = mwRandomUnitPoint(dsfmtState);
        tmp = mw_sqrv(&vec);
        rsq = showRealValue(&tmp);         /* compute radius squared */
    }
    while (rsq > 1.0);              /* reject if outside sphere */

    rsc = rad / mw_sqrt_0(rsq);       /* compute scaling factor */

    vec.x = mw_mul_s(&vec.x, rsc);
    vec.y = mw_mul_s(&vec.y, rsc);
    vec.z = mw_mul_s(&vec.z, rsc);    /* rescale to radius given */

    return vec;
}

static real_0 nfwRandomR(dsfmt_t* dsfmtState, real_0 startradius, real_0 endradius)
{
    real_0 rnd;

    /* returns [0, 1) */
    rnd = (real_0) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    return (endradius - startradius) * rnd + startradius;
}

static real_0 nfwSelectFromG(dsfmt_t* dsfmtState)
{
    real_0 x, y;

    do                      /* select from fn g(x) */
    {
        x = mwXrandom(dsfmtState, 0.0, 1.0);      /* for x in range 0:1 */
        y = mwXrandom(dsfmtState, 0.0, 0.1);      /* max of g(x) is 0.092 */
    }   /* using von Neumann tech */
    while (y > sqr_0(x) * mw_pow_0(1.0 - sqr_0(x), 3.5));

    return x;
}

static real_0 nfwCalculateV(real_0 r, real_0 rho_0, real_0 R_S)
{
    real_0 v;
    real_0 mass = nfwMassInsideRadius(r, rho_0, R_S);
    v = mw_sqrt_0( /*G!!!*/ mass / r);

    return v;
}

static mwvector nfwBodyPosition(dsfmt_t* dsfmtState, mwvector* rshift, real_0 rsc, real_0 r)
{
    mwvector pos;

    pos = nfwPickShell(dsfmtState, rsc * r);  /* pick scaled position */
    pos = mw_addv(&pos, rshift);               /* move the position */

    return pos;
}

static mwvector nfwBodyVelocity(dsfmt_t* dsfmtState, mwvector* vshift, real_0 r, real_0 rho_0, real_0 R_S)
{
    mwvector vel;
    real_0 v;

    v = nfwCalculateV(r, rho_0, R_S);
    vel = nfwPickShell(dsfmtState, v);   /* pick scaled velocity */
    vel = mw_addv(&vel, vshift);             /* move the velocity */

    return vel;
}

/* generatenfw: generate nfw model initial conditions
 * Extremely hacky. If you actually want to use this
 * talk to Colin Rice before you do anything. Seriously.
 */
static int nbGenerateNFWCore(lua_State* luaSt,

                             dsfmt_t* prng,
                             unsigned int nbody,
                             real_0 mass,

                             mwbool ignore,

                             mwvector* rShift,
                             mwvector* vShift,
                             real_0 rho_0,
                             real_0 R_S)
{
    unsigned int i;
    int table;
    Body b;
    real_0 r;
    real_0 totalMass = 0.0;
    real_0 radius = 0.0;
    real_0 massEpsilon = mass / nbody;    /* The amount of mass we increase for
                                          each particle */

    memset(&b, 0, sizeof(b));

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mw_real_const(mass / nbody);    /* Mass per particle */


    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);


    /* Start with half an epsilon */
    totalMass = 0.5 * massEpsilon;

    for (i = 0; i < nbody; ++i)
    {
        real_0 endradius = nfwNextRadius(radius, totalMass + massEpsilon, rho_0, R_S);

        do
        {
            r = nfwRandomR(prng, radius, endradius);
        }
        while (isinf(r));

        radius = endradius;
        totalMass += massEpsilon;

        b.bodynode.pos = nfwBodyPosition(prng, rShift, 1, r);
        b.vel = nfwBodyVelocity(prng, vShift, r, rho_0, R_S);
        assert(nbPositionValid(&b.bodynode.pos));

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

int nbGenerateNFW(lua_State* luaSt)
{
    static dsfmt_t* prng;
    static const mwvector* position = NULL;
    static const mwvector* velocity = NULL;
    static mwbool ignore;
    static real_0 mass = 0.0, nbodyf = 0.0, rho_0 = 0.0, R_S = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "nbody",        LUA_TNUMBER,   NULL,          TRUE,  &nbodyf   },
            { "mass",         LUA_TNUMBER,   NULL,          TRUE,  &mass     },
            { "rho_0",        LUA_TNUMBER,   NULL,          TRUE,  &rho_0    },
            { "scaledRadius", LUA_TNUMBER,   NULL,          TRUE,  &R_S      },
            { "position",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &position },
            { "velocity",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &velocity },
            { "ignore",       LUA_TBOOLEAN,  NULL,          FALSE, &ignore   },
            { "prng",         LUA_TUSERDATA, DSFMT_TYPE,    TRUE,  &prng     },
            END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 arguments");

    handleNamedArgumentTable(luaSt, argTable, 1);

    return nbGenerateNFWCore(luaSt, prng, (unsigned int) nbodyf, mass, ignore,
                                 position, velocity, rho_0, R_S);
}

void registerGenerateNFW(lua_State* luaSt)
{
    lua_register(luaSt, "generateNFW", nbGenerateNFW);
}

