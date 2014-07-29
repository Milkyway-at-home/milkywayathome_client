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

/* pickshell: pick a random point on a sphere of specified radius. 
*
*   Changed this section to be in compliance with NEMO's initialization technique.
*   Instead of assigning it a point on sphere, assigns angles. Allows for non-circular orbits.
*
*/
static inline mwvector pickShell(dsfmt_t* dsfmtState, real rad)
{
    real rsq, rsc;
    mwvector vec;
    real phi, theta;

    /*defining some angles*/
    theta = mw_acos( mwXrandom(dsfmtState, -1.0, 1.0) );
    phi =   mwXrandom( dsfmtState, 0.0, 2.0*M_PI );

    /*this is standard formula for x,y,z components in spherical*/
    X(vec) = rad*sin( theta )*cos( phi );    /*x component*/
    Y(vec) = rad*sin( theta )*sin( phi );    /*y component*/
    Z(vec) = rad*cos( theta );               /*z component*/

    rsq = mw_sqrv(vec);             /* compute radius squared */
    rsc = rad / mw_sqrt(rsq);       /* compute scaling factor */
    mw_incmulvs(vec, rsc);          /* rescale to radius given */

    return vec;
}

/*
*   This code snippet is from nbody_plummer.c.
*   This is a probability distribution which is used in NEMO.
*   The parameter that is returned is a fraction that is used to sample the velocity.
*   See Aarseth et al. (1974), eq. (A4,5).
*/
static inline real plummerSelectFromG(dsfmt_t* dsfmtState)
{

    real x, y;

    do                      /* select from fn g(x) */
    {
        x = mwXrandom(dsfmtState, 0.0, 1.0);      /* for x in range 0:1 */
        y = mwXrandom(dsfmtState, 0.0, 0.1);      /* max of g(x) is 0.092 */
    }   /* using von Neumann tech */
    while (y > x*x * mw_pow(1.0 - x*x, 3.5));

    return x;
}


static inline real isotropicRandomR(dsfmt_t* dsfmtState, real scaleRad1, real scaleRad2,
				    real Mass1, real Mass2)
{
  // Rejection sampling radius generation
  real scaleRad1Cube = cube(scaleRad1);
  real scaleRad2Cube = cube(scaleRad2);
  real RHO_MAX = 3.0/(4.0*M_PI) * (Mass1/scaleRad1Cube + Mass2/scaleRad2Cube);
  mwbool GOOD_RADIUS = 0;
  // Arbitrarily define sample range to be [0, 5(a1 + a2)]

  real r;
  real u, val;

  while (GOOD_RADIUS != 1)
    {
      r = mwXrandom(dsfmtState,0.0, 3.0 * (scaleRad1 + scaleRad2));
      u = (real)mwXrandom(dsfmtState,0.0,1.0);

      val = 3.0/(4.0 *M_PI)*(Mass1/scaleRad1Cube * mw_pow(1.0 + sqr(r)/ sqr(scaleRad1),-2.5) +
			     Mass2/scaleRad2Cube * mw_pow(1.0 + sqr(r)/sqr(scaleRad2),-2.5));

      if (val/RHO_MAX > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  return r;
}

/*
*   Added the x value. This is a parameter that is used to sample the velocity.
*   Though, the distribution function was originally used for one plummer sphere.
*   We have two here. Though I think the distribution function can still be used.
*/
static inline real isotropicRandomV(dsfmt_t* dsfmtState,real r, real scaleRad1, real scaleRad2, real Mass1, real Mass2)
{

  real GMsolar = mw_pow(1.327,20.0); //SI UNITS m^3 / sec^2
  scaleRad1 *= 3.086e19; //meters
  scaleRad2 *= 3.086e19;  
  r *= 3.086e19;
  real val;
  real x;
  x= plummerSelectFromG(dsfmtState);

  val = x * M_SQRT2* mw_sqrt(GMsolar *Mass1/mw_sqrt(sqr(r) + sqr(scaleRad1))
			       + GMsolar * Mass2/mw_sqrt(sqr(r) + sqr(scaleRad2)));

  return val; //km/s
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

    v = isotropicRandomV(dsfmtState,r,scaleRad1,scaleRad2,Mass1,Mass2);
    vel = pickShell(dsfmtState, vsc * v);   /* pick scaled velocity */
    mw_incaddv(vel, vshift);                /* move the velocity */

    return vel;
}

/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Heggie,
 * etc).  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */

/*
* I have chosen not to use Henon units for the reason that the distribution
* function is not easily decomposed into a function of the system mass and
* Newton's gravitational constant.  Instead, SI units are used and then
* converted to km/s for the velocities at the end. Properly, one should use
* CGS units as this would conform to the standards of theoretical astrophysics,
* but I am lazy and leave this as an exercise to a future graduate student :)))
*
* Lots of love and rainbows,
*
* Jake B.
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
    memset(&b, 0, sizeof(b));

    velScale = 1000;// Conversion from km/s

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


