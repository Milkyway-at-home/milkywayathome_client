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
#include "nbody_plummer.h"

/* pickshell: pick a random point on a sphere of specified radius. */
static inline mwvector pickShell(dsfmt_t* dsfmtState, real* rad)
{
    real rsq, rsc, tmp;
    mwvector vec;

    do                      /* pick point in NDIM-space */
    {
        vec = mwRandomUnitPoint(dsfmtState);
        rsq = mw_sqrv(&vec);         /* compute radius squared */
    }
    while (showRealValue(&rsq) > 1.0);              /* reject if outside sphere */

    tmp = mw_sqrt(&rsq);
    rsc = mw_div(rad, &tmp);       /* compute scaling factor */

    vec.x = mw_mul(&vec.x, &rsc);          /* rescale to radius given */
    vec.y = mw_mul(&vec.y, &rsc);
    vec.z = mw_mul(&vec.z, &rsc);

    return vec;
}

static inline real_0 plummerRandomR(dsfmt_t* dsfmtState)
{
    real_0 rnd;

    /* returns [0, 1) */
    rnd = (real_0) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    return 1.0 / mw_sqrt_0(mw_pow_0(rnd, -2.0 / 3.0) - 1.0);
}

static inline real_0 plummerSelectFromG(dsfmt_t* dsfmtState)
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

static inline real_0 plummerRandomV(dsfmt_t* dsfmtState, real_0 r)
{
    real_0 x, v;

    x = plummerSelectFromG(dsfmtState);
    v = M_SQRT2 * x / mw_sqrt_0(mw_sqrt_0(1.0 + sqr_0(r)));   /* find v in struct units */

    return v;
}

static inline mwvector plummerBodyPosition(dsfmt_t* dsfmtState, mwvector* rshift, real* rsc, real_0 r)
{
    mwvector pos;
    real tmp = mw_mul_s(rsc, r);

    pos = pickShell(dsfmtState, &tmp);  /* pick scaled position */
    pos = mw_addv(&pos, rshift);               /* move the position */

    return pos;
}

static inline mwvector plummerBodyVelocity(dsfmt_t* dsfmtState, mwvector* vshift, real* vsc, real_0 r)
{
    mwvector vel;
    real_0 v;

    v = plummerRandomV(dsfmtState, r);
    real tmp = mw_mul_s(vsc, v);
    vel = pickShell(dsfmtState, &tmp);   /* pick scaled velocity */
    vel = mw_addv(&vel, vshift);                /* move the velocity */
    
    return vel;
}

/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Hegge,
 * etc).  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */
static int nbGeneratePlummerCore(lua_State* luaSt,

                                 dsfmt_t* prng,
                                 unsigned int nbody,
                                 real* mass,

                                 mwbool ignore,

                                 mwvector* rShift,
                                 mwvector* vShift,
                                 real* radiusScale)
{
    unsigned int i;
    int table;
    Body b;
    real_0 r;
    real velScale, tmp;

    memset(&b, 0, sizeof(b));

    tmp = mw_div(&mass, &radiusScale);
    velScale = mw_sqrt(&tmp);     /* and recip. speed scale */

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mw_mul_s(mass, inv_0((real_0) nbody));    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < nbody; ++i)
    {
        do
        {
            r = plummerRandomR(prng);
            /* FIXME: We should avoid the divide by 0.0 by multiplying
             * the original random number by 0.9999.. but I'm too lazy
             * to change the tests. Same with other models */
        }
        while (isinf(r));
        
        b.bodynode.id = i + 1;
        b.bodynode.pos = plummerBodyPosition(prng, rShift, radiusScale, r);
        //mw_printf("POS = [ %.15f, %.15f, %.15f ]\n", showRealValue(&b.bodynode.pos.x), showRealValue(&b.bodynode.pos.y), showRealValue(&b.bodynode.pos.z));
        b.vel = plummerBodyVelocity(prng, vShift, &velScale, r);
        //mw_printf("VEL = [ %.15f, %.15f, %.15f ]\n", showRealValue(&b.vel.x), showRealValue(&b.vel.y), showRealValue(&b.vel.z));

        assert(nbPositionValid(&b.bodynode.pos));

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

int nbGeneratePlummer(lua_State* luaSt)
{
    static dsfmt_t* prng;
    static const mwvector* position = NULL;
    static const mwvector* velocity = NULL;
    static mwbool ignore;
    static real_0 mass = 0.0;
    static real_0 nbodyf = 0.0;
    static real_0 radiusScale = 1.0;

    static const MWNamedArg argTable[] =
        {
            { "nbody",        LUA_TNUMBER,   NULL,          TRUE,  &nbodyf      },
            { "prng",         LUA_TUSERDATA, DSFMT_TYPE,    TRUE,  &prng        },
            { "position",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &position    },
            { "velocity",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &velocity    },
            { "mass",         LUA_TNUMBER,   NULL,          TRUE,  &mass        },
            { "scaleRadius",  LUA_TNUMBER,   NULL,          TRUE,  &radiusScale },
            { "ignore",       LUA_TBOOLEAN,  NULL,          FALSE, &ignore      },
            END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 arguments");

    handleNamedArgumentTable(luaSt, argTable, 1);

    real real_mass = mw_real_var(mass, BARYON_MASS_POS);
    real real_rad = mw_real_var(radiusScale, BARYON_RADIUS_POS);

    return nbGeneratePlummerCore(luaSt, prng, (unsigned int) nbodyf, &real_mass, ignore,
                                 position, velocity, &real_rad);
}

void registerGeneratePlummer(lua_State* luaSt)
{
    lua_register(luaSt, "generatePlummer", nbGeneratePlummer);
}

