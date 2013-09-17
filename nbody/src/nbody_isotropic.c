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
#include "milkyway_math.h"
#include "milkyway_lua.h"
#include "nbody_lua_types.h"
#include "nbody_isotropic.h"

/* pickshell: pick a random point on a sphere of specified radius. */
static inline mwvector pickShell(dsfmt_t* dsfmtState, real rad)
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

static inline real isotropicRandomR(dsfmt_t* dsfmtState, real scaleRad1, real scaleRad2,
				    real Mass1, real Mass2)
{
  // Rejection sampling radius generation
  real RHO_MAX = 3/(4*M_PI) * (Mass1/(mw_pow(scaleRad1,3)) + Mass2/(mw_pow(scaleRad2,3)));
  mwbool GOOD_RADIUS = 0;
  // Arbitrarily define sample range to be [0, 5(a1 + a2)

  real r;

  while (GOOD_RADIUS != 1)
    {
      r = mwXrandom(dsfmtState,0.0, 3*(scaleRad1 + scaleRad2));
      real u = (real)mwXrandom(dsfmtState,0.0,1.0);

      real val = 3/(4*M_PI)*(Mass1/(mw_pow(scaleRad1,3)) *mw_pow(1 + mw_pow(r,2)/mw_pow(scaleRad1,2),-5/2) + 
			     Mass2/(mw_pow(scaleRad2,3))*mw_pow(1 + mw_pow(r,2)/mw_pow(scaleRad2,2),-5/2));
      
      if (val/RHO_MAX > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  return r;
}



static inline real isotropicRandomV(real r, real scaleRad1, real scaleRad2,                                                                                                             				    real Mass1, real Mass2)  
{


  real val;
  val = mw_sqrt(Mass1/mw_sqrt(mw_pow(r,2) + mw_pow(scaleRad1,2)) 
		+ Mass2/mw_sqrt(mw_pow(r,2) + mw_pow(scaleRad2,2)));

  return val;
}

static inline mwvector isotropicBodyPosition(dsfmt_t* dsfmtState, mwvector rshift,  real r)
{
    mwvector pos;

    pos = pickShell(dsfmtState,  r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */

    return pos;
}

static inline mwvector isotropicBodyVelocity(dsfmt_t* dsfmtState,real r, mwvector vshift, real vsc,  real scaleRad1, real scaleRad2,
					     real Mass1, real Mass2)
{
    mwvector vel;
    real v;

    v = isotropicRandomV(r,scaleRad1,scaleRad2,Mass1,Mass2);
    vel = pickShell(dsfmtState, vsc * v);   /* pick scaled velocity */
    mw_incaddv(vel, vshift);                /* move the velocity */

    return vel;
}

/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Heggie,
 * etc).  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */
static int nbGenerateIsotropicCore(lua_State* luaSt,

				   dsfmt_t* prng,
				   unsigned int nbody,
				   real mass1,
				   real mass2,

				   mwbool ignore,
				   
				   mwvector rShift,
				   mwvector vShift,
				   real radiusScale1,
				   real radiusScale2)
{
    unsigned int i;
    int table;
    Body b;
    real r, velScale;

    real mass = mass1 + mass2;
    real radiusScale = mw_sqrt(mw_pow(radiusScale1,2) + mw_pow(radiusScale2,2));
    memset(&b, 0, sizeof(b));

    velScale =  mw_sqrt(mass / radiusScale);     /* and recip. speed scale */

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);
    
    for (i = 0; i < nbody; ++i)
    {
        do
        {
	  r = isotropicRandomR(prng, radiusScale1, radiusScale2, mass1, mass2);
	    /* FIXME: We should avoid the divide by 0.0 by multiplying
             * the original random number by 0.9999.. but I'm too lazy
             * to change the tests. Same with other models */
        }
        while (isinf(r));

        b.bodynode.pos = isotropicBodyPosition(prng, rShift, r);

        b.vel = isotropicBodyVelocity(prng, r, vShift, velScale, radiusScale1, radiusScale2, mass1, mass2);
	
        assert(nbPositionValid(b.bodynode.pos));

        pushBody(luaSt, &b);
	//	printf("Body %d is pushed. \n",i);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

int nbGenerateIsotropic(lua_State* luaSt)
{
    static dsfmt_t* prng;
    static const mwvector* position = NULL;
    static const mwvector* velocity = NULL;
    static mwbool ignore;
    static real mass1 = 0.0, nbodyf = 0.0, radiusScale1 = 0.0;
    static real mass2 = 0.0, radiusScale2 = 0.0;

    static const MWNamedArg argTable[] =
        {
	  { "nbody",        LUA_TNUMBER,   NULL,          TRUE,  &nbodyf      },
	  { "mass1",        LUA_TNUMBER,   NULL,          TRUE,  &mass1       },
          { "mass2",        LUA_TNUMBER,   NULL,          TRUE,  &mass2       },
	  { "scaleRadius1", LUA_TNUMBER,   NULL,          TRUE,  &radiusScale1},
	  { "scaleRadius2", LUA_TNUMBER,   NULL,          TRUE,  &radiusScale2},
	  { "position",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &position    },
	  { "velocity",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &velocity    },
	  { "ignore",       LUA_TBOOLEAN,  NULL,          FALSE, &ignore      },
	  { "prng",         LUA_TUSERDATA, DSFMT_TYPE,    TRUE,  &prng        },
	  END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 arguments");

    handleNamedArgumentTable(luaSt, argTable, 1);

    return nbGenerateIsotropicCore(luaSt, prng, (unsigned int) nbodyf, mass1, mass2, ignore,
                                 *position, *velocity, radiusScale1, radiusScale2);
}

void registerGenerateIsotropic(lua_State* luaSt)
{
    lua_register(luaSt, "generateIsotropic", nbGenerateIsotropic);
}

