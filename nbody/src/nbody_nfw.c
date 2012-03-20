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

real nfwMassInsideRadius(real radius, real rho_0, real R_S){
    //Returns that mass inside a certain radius

    //Leading Constant
    real mass = 4.0 * M_PI * rho_0 * cube(radius);
    //Integration terms
    mass *= ( mw_log((R_S + radius) / R_S) - radius / (R_S + radius));
    return mass;

}

real nfwNextRadius(real start_radius, real goal_mass, real rho_0, real R_S){
    // This is a scary function which returns the next radius limit
    // Do not use this code on client computers. It is slow inneficient
    // and well just dirty.

    for (real test_radius = start_radius;; test_radius += (real)0.0001){
        if (nfwMassInsideRadius(test_radius, rho_0, R_S) >= goal_mass)
        {
            return test_radius;
        }
    }
}

/* nfwPickShell: pick a random point on a sphere of specified radius. */
mwvector nfwPickShell(dsfmt_t* dsfmtState, real rad)
{
    real rsq, rsc;
    mwvector vec;

    do                      /* pick point in NDIM-space */
    {
        vec = mwRandomUnitPoint(dsfmtState);
        rsq = mw_sqrv(vec);         /* compute radius squared */
    }
    while (rsq > 1.0);              /* reject if outside sphere */

    rsc = rad / mw_sqrt(rsq);       /* compute scaling factor */
    mw_incmulvs(vec, rsc);          /* rescale to radius given */

    return vec;
}

static inline real nfwRandomR(dsfmt_t* dsfmtState, real startradius, real endradius)
{
    real rnd;

    /* returns [0, 1) */
    rnd = (real) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    return (endradius - startradius) * rnd + startradius;
}

static inline real nfwSelectFromG(dsfmt_t* dsfmtState)
{
    real x, y;

    do                      /* select from fn g(x) */
    {
        x = mwXrandom(dsfmtState, 0.0, 1.0);      /* for x in range 0:1 */
        y = mwXrandom(dsfmtState, 0.0, 0.1);      /* max of g(x) is 0.092 */
    }   /* using von Neumann tech */
    while (y > sqr(x) * mw_pow(1.0 - sqr(x), 3.5));

    return x;
}

static inline real nfwCalculateV(real r, real rho_0, real R_S)
{
    real v;
    real mass = nfwMassInsideRadius(r, rho_0, R_S);
    v = mw_sqrt( /*G!!!*/ mass / r);

    return v;
}

static inline mwvector nfwBodyPosition(dsfmt_t* dsfmtState, mwvector rshift, real rsc, real r)
{
    mwvector pos;

    pos = nfwPickShell(dsfmtState, rsc * r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */

    return pos;
}

static inline mwvector nfwBodyVelocity(dsfmt_t* dsfmtState, mwvector vshift, real r, real rho_0, real R_S)
{
    mwvector vel;
    real v;

    v = nfwCalculateV(r, rho_0, R_S);
    vel = nfwPickShell(dsfmtState, v);   /* pick scaled velocity */
    mw_incaddv(vel, vshift);              /* move the velocity */

    return vel;
}

/* generatenfw: generate nfw model initial conditions
 * Extremely hacky. If you actually want to use this
 * talk to Colin Rice before you do anything. Seriously.
 */
static int nbGenerateNFWCore(lua_State* luaSt,

                                 dsfmt_t* prng,
                                 unsigned int nbody,
                                 real mass,

                                 mwbool ignore,

                                 mwvector rShift,
                                 mwvector vShift,
                                 real rho_0,
                                 real R_S)
{
    unsigned int i;
    int table;
    Body b;
    real r;

    memset(&b, 0, sizeof(b));

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */
    real massEpsilon = mass / nbody;    /* The amount of mass we increase for
                                        each particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);

    real radius = 0;

    //Start with half an epsilon
    real totalMass = massEpsilon / (real)2;

    for (i = 0; i < nbody; ++i)
    {
        real endradius = nfwNextRadius(radius, totalMass + massEpsilon, rho_0, R_S);

        r = nfwRandomR(prng, radius, endradius);

        radius = endradius;
        totalMass += massEpsilon;

        b.bodynode.pos = nfwBodyPosition(prng, rShift, 1, r);
        b.vel = nfwBodyVelocity(prng, vShift, r, rho_0, R_S);

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
    static real mass = 0.0, nbodyf = 0.0, rho_0 = 0.0, R_S = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "nbody",        LUA_TNUMBER,   NULL,          TRUE,  &nbodyf      },
            { "mass",         LUA_TNUMBER,   NULL,          TRUE,  &mass        },
            { "rho_0",        LUA_TNUMBER,   NULL,          TRUE,     &rho_0
             },
            { "scaledRadius",          LUA_TNUMBER,   NULL,          TRUE,     &R_S
             },
            { "position",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &position    },
            { "velocity",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &velocity    },
            { "ignore",       LUA_TBOOLEAN,  NULL,          FALSE, &ignore      },
            { "prng",         LUA_TUSERDATA, DSFMT_TYPE,    TRUE,  &prng        },
            END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 arguments");

    handleNamedArgumentTable(luaSt, argTable, 1);

    return nbGenerateNFWCore(luaSt, prng, (unsigned int) nbodyf, mass, ignore,
                                 *position, *velocity, rho_0, R_S);
}

void registerGenerateNFW(lua_State* luaSt)
{
    lua_register(luaSt, "generateNFW", nbGenerateNFW);
}

