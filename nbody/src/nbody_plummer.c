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
#include "nbody_util.h"
#include "milkyway_util.h"
#include "milkyway_lua_marshal.h"
#include "nbody_lua_types.h"
#include "nbody_plummer.h"

/* pickshell: pick a random point on a sphere of specified radius. */
static inline mwvector pickShell(dsfmt_t* dsfmtState, real rad)
{
    real rsq, rsc;
    mwvector vec;

    do                      /* pick point in NDIM-space */
    {
        vec = mwRandomVector(dsfmtState);
        rsq = mw_sqrv(vec);         /* compute radius squared */
    }
    while (rsq > 1.0);              /* reject if outside sphere */

    rsc = rad / mw_sqrt(rsq);       /* compute scaling factor */
    mw_incmulvs(vec, rsc);          /* rescale to radius given */

    return vec;
}

static void printPlummer(mwvector rshift, mwvector vshift)
{
    warn("<plummer_r> %.14g %.14g %.14g </plummer_r>\n"
         "<plummer_v> %.14g %.14g %.14g </plummer_v>\n",
         X(rshift), Y(rshift), Z(rshift),
         X(vshift), Y(vshift), Z(vshift));
}

static inline real plummerRandomR(dsfmt_t* dsfmtState)
{
    real rnd;

    /* returns [0, 1) */
    rnd = (real) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    return 1.0 / mw_sqrt(mw_pow(rnd, -2.0 / 3.0) - 1.0);
}

static inline real plummerSelectFromG(dsfmt_t* dsfmtState)
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

static inline real plummerRandomV(dsfmt_t* dsfmtState, real r)
{
    real x, v;

    x = plummerSelectFromG(dsfmtState);
    v = M_SQRT2 * x / mw_sqrt(mw_sqrt(1.0 + sqr(r)));   /* find v in struct units */

    return v;
}

static inline mwvector plummerBodyPosition(dsfmt_t* dsfmtState, mwvector rshift, real rsc, real r)
{
    mwvector pos;

    pos = pickShell(dsfmtState, rsc * r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */

    return pos;
}

static inline mwvector plummerBodyVelocity(dsfmt_t* dsfmtState, mwvector vshift, real vsc, real r)
{
    mwvector vel;
    real v;

    v = plummerRandomV(dsfmtState, r);
    vel = pickShell(dsfmtState, vsc * v);   /* pick scaled velocity */
    mw_incaddv(vel, vshift);                /* move the velocity */

    return vel;
}

/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Hegge,
 * etc).  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */
static int createPlummerSphereTable(lua_State* luaSt,
                                    dsfmt_t* prng,
                                    unsigned int nbody,
                                    real mass,

                                    mwbool ignore,

                                    mwvector rShift,
                                    mwvector vShift,
                                    real radiusScale)
{
    unsigned int i;
    int table;
    body b;
    real r, velScale;

    velScale = mw_sqrt(mass / radiusScale);     /* and recip. speed scale */

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);
    for (i = 0; i < nbody; ++i)
    {
        r = plummerRandomR(prng);

        b.bodynode.pos = plummerBodyPosition(prng, rShift, radiusScale, r);
        b.vel = plummerBodyVelocity(prng, vShift, velScale, r);

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

int generatePlummer(lua_State* luaSt)
{
    static dsfmt_t* prng;
    static const InitialConditions* ic;
    static mwbool ignore;
    static real mass;
    static int nbody;
    static real radiusScale = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "prng",              LUA_TUSERDATA, DSFMT_TYPE,              TRUE,  &prng        },
            { "initialConditions", LUA_TUSERDATA, INITIAL_CONDITIONS_TYPE, TRUE,  &ic          },
            { "scaleRadius",       LUA_TNUMBER,   NULL,                    TRUE,  &radiusScale },
            { "mass",              LUA_TNUMBER,   NULL,                    TRUE,  &mass        },
            { "ignore",            LUA_TBOOLEAN,  NULL,                    FALSE, &ignore      },
            END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 2)
        return luaL_argerror(luaSt, 1, "Expected 2 arguments");

    nbody = luaL_checkinteger(luaSt, 1);
    handleNamedArgumentTable(luaSt, argTable, 2);

    return createPlummerSphereTable(luaSt, prng, nbody, mass, ignore,
                                    ic->position, ic->velocity, radiusScale);
}

void registerGeneratePlummer(lua_State* luaSt)
{
    lua_register(luaSt, "generatePlummer", generatePlummer);
}



