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
#include "lua_type_marshal.h"
#include "nbody_lua_types.h"

static mwvector randomVec(dsfmt_t* dsfmtState)
{
    /* pick from unit cube */
    mwvector vec;

    X(vec) = mwUnitRandom(dsfmtState);
    Y(vec) = mwUnitRandom(dsfmtState);
    Z(vec) = mwUnitRandom(dsfmtState);
    W(vec) = 0.0;

    return vec;
}

/* pickshell: pick a random point on a sphere of specified radius. */
static inline mwvector pickShell(dsfmt_t* dsfmtState, real rad)
{
    real rsq, rsc;
    mwvector vec;

    do                      /* pick point in NDIM-space */
    {
        vec = randomVec(dsfmtState);
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
    real rnd, r;

    /* returns [0, 1) */
    rnd = (real) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    r = 1.0 / mw_sqrt(mw_pow(rnd, -2.0 / 3.0) - 1.0);

    return r;
}

static inline real plummerSelectFromG(dsfmt_t* dsfmtState)
{
    real x, y;

    do                      /* select from fn g(x) */
    {
        x = mwXrandom(dsfmtState, 0.0, 1.0);      /* for x in range 0:1 */
        y = mwXrandom(dsfmtState, 0.0, 0.1);      /* max of g(x) is 0.092 */
    }   /* using von Neumann tech */
    while (y > -cube(x - 1.0) * sqr(x) * cube(x + 1.0) * mw_sqrt(1.0 - sqr(x)));

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
                                    real radiusScale,
                                    real velScale)
{
    unsigned int i;
    int table;
    body b;
    real r;

    Type(&b) = BODY(ignore);    /* Same for all in the model */
    Mass(&b) = mass / nbody;    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);
    for (i = 0; i < nbody; ++i)
    {
        r = plummerRandomR(prng);

        Pos(&b) = plummerBodyPosition(prng, rShift, radiusScale, r);
        Vel(&b) = plummerBodyVelocity(prng, vShift, velScale, r);

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

/* DSFMT -> InitialConditions -> Int (nbody) -> Real (mass) -> Bool (ignore final) -> [UserData] -> [Body] */
int generatePlummer(lua_State* luaSt)
{
    int nArgs, userDataIndex;
    dsfmt_t* prng;
    const InitialConditions* ic;
    mwbool ignore;
    real mass, velScale;
    int nbody;
    static real radiusScale = 0.0;

    static const MWNamedArg udTable[] =
        {
            { "scaleRadius", LUA_TNUMBER, NULL, TRUE, &radiusScale },  /* length scale factor */
            END_MW_NAMED_ARG
        };

    nArgs = lua_gettop(luaSt);

    switch (nArgs)
    {
        case 6:
            prng = checkDSFMT(luaSt, 1);
            ic = checkInitialConditions(luaSt, 2);
            nbody = luaL_checkinteger(luaSt, 3);
            mass = luaL_checknumber(luaSt, 4);
            ignore = mw_lua_checkboolean(luaSt, 5);
            userDataIndex = mw_lua_checktable(luaSt, 6);

            handleNamedArgumentTable(luaSt, udTable, userDataIndex);
            break;

        default:
            return luaL_argerror(luaSt, 1, "Expected 6 arguments");
    }

    velScale = mw_sqrt(mass / radiusScale);     /* and recip. speed scale */

    return createPlummerSphereTable(luaSt, prng, nbody, mass, ignore,
                                    ic->position, ic->velocity, radiusScale, velScale);
}

void registerGeneratePlummer(lua_State* luaSt)
{
    lua_pushcfunction(luaSt, generatePlummer);
    lua_setglobal(luaSt, "generatePlummer");
}



