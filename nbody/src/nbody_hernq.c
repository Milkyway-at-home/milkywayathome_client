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
#include "nbody_hernq.h"

//FIXME: THIS CODE DOES NOT PROPAGATE DERIVAIVE INFORMATION! MUST BE REWORKED BEFORE RUNNING WITH AUTODIFF!

static real_0 hernqMassInsideRadius(real_0 radius, real_0 radius_scale, real_0 a, real_0 mass)
{
    //Returns that mass inside a certain radius
    radius = radius / radius_scale;

    //Integration terms
    mass = mass * radius * radius / ((a + radius) * (a + radius)) ;
    return mass;

}

static real_0 hernqNextRadius(real_0 startRadius, real_0 goalMass, real_0 radius, real_0 a, real_0 mass)
{
    // This is a scary function which returns the next radius limit
    // Do not use this code on client computers. It is slow inneficient
    // and well just dirty.

    real_0 test_radius;

    for (test_radius = startRadius;; test_radius += (real_0)0.0001)
    {
        if (hernqMassInsideRadius(test_radius, radius, a, mass) >= goalMass)
        {
            return test_radius;
        }
    }
}

/* hernqPickShell: pick a random point on a sphere of specified radius. */
static mwvector hernqPickShell(dsfmt_t* dsfmtState, real_0 rad)
{
    real rsq, rsc, tmp;
    mwvector vec;

    do                      /* pick point in NDIM-space */
    {
        vec = mwRandomUnitPoint(dsfmtState);
        rsq = mw_sqrv(&vec);         /* compute radius squared */
    }
    while (showRealValue(&rsq) > 1.0);              /* reject if outside sphere */

    tmp = minushalf(&rsq);
    rsc = mw_mul_s(&tmp, rad);       /* compute scaling factor */
    mw_incmulvs(&vec, &rsc);          /* rescale to radius given */

    return vec;
}

static real_0 hernqRandomR(dsfmt_t* dsfmtState, real_0 startradius, real_0 endradius)
{
    real_0 rnd;

    /* returns [0, 1) */
    rnd = (real_0) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    return (endradius - startradius) * rnd + startradius;
}

static real_0 hernqSelectFromG(dsfmt_t* dsfmtState)
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

static real_0 hernqCalculateV(real_0 r, real_0 radius, real_0 a, real_0 mass)
{
    real_0 v;
    mass = hernqMassInsideRadius(r, radius, a, mass);
    v = mw_sqrt_0( /*G!!!*/ mass / r);

    return v;
}

static mwvector hernqBodyPosition(dsfmt_t* dsfmtState, mwvector* rshift, real_0 rsc, real_0 r)
{
    mwvector pos;

    pos = hernqPickShell(dsfmtState, rsc * r);  /* pick scaled position */
    mw_incaddv(&pos, rshift);               /* move the position */

    return pos;
}

static mwvector hernqBodyVelocity(dsfmt_t* dsfmtState, mwvector* vshift, real_0 r, real_0 radius, real_0 a, real_0 mass)
{
    mwvector vel;
    real_0 v;

    v = hernqCalculateV(r, radius, a, mass);
    vel = hernqPickShell(dsfmtState, v);   /* pick scaled velocity */
    mw_incaddv(&vel, vshift);              /* move the velocity */

    return vel;
}

/* generateHernq: generate hernquist model initial conditions
 * Extremely hacky. If you actually want to use this
 * talk to Colin Rice before you do anything. Seriously.
 */
static int nbGenerateHernqCore(lua_State* luaSt,

                                 dsfmt_t* prng,
                                 unsigned int nbody,
                                 real_0 mass,

                                 mwbool ignore,

                                 mwvector* rShift,
                                 mwvector* vShift,
                                 real_0 radius_scale,
                                 real_0 a)
{
    unsigned int i;
    int table;
    Body b;
    real_0 r;
    real_0 radius = 0.0;
    real_0 massEpsilon = mass / nbody; /* The amount of mass we increase for
                                        each particle */

    /* Start with half an epsilon */
    real_0 totalMass = 0.5 * massEpsilon;

    memset(&b, 0, sizeof(b));

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mw_real_const(mass / nbody);    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < nbody; ++i)
    {
        real_0 endradius = hernqNextRadius(radius, totalMass + massEpsilon, radius, a, mass);

        do
        {
            r = hernqRandomR(prng, radius, endradius);
        }
        while (isinf(r));

        radius = endradius;
        totalMass += massEpsilon;

        b.bodynode.pos = hernqBodyPosition(prng, rShift, 1, r);
        b.vel = hernqBodyVelocity(prng, vShift, r, radius_scale, a, mass);
        assert(nbPositionValid(&b.bodynode.pos));

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

int nbGenerateHernq(lua_State* luaSt)
{
    static dsfmt_t* prng;
    static const mwvector* position = NULL;
    static const mwvector* velocity = NULL;
    static mwbool ignore;
    static real_0 mass = 0.0, nbodyf = 0.0, radius = 0.0, a = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "nbody",    LUA_TNUMBER,   NULL,          TRUE,  &nbodyf   },
            { "mass",     LUA_TNUMBER,   NULL,          TRUE,  &mass     },
            { "radius",   LUA_TNUMBER,   NULL,          TRUE,  &radius   },
            { "a",        LUA_TNUMBER,   NULL,          TRUE,  &a        },
            { "position", LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &position },
            { "velocity", LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &velocity },
            { "ignore",   LUA_TBOOLEAN,  NULL,          FALSE, &ignore   },
            { "prng",     LUA_TUSERDATA, DSFMT_TYPE,    TRUE,  &prng     },
            END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 arguments");

    handleNamedArgumentTable(luaSt, argTable, 1);

    return nbGenerateHernqCore(luaSt, prng, (unsigned int) nbodyf, mass, ignore,
                                 position, velocity, radius, a);
}

void registerGenerateHernq(lua_State* luaSt)
{
    lua_register(luaSt, "generateHernq", nbGenerateHernq);
}

