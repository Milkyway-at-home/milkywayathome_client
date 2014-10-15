/*THIS FILE IS IN DEVELOPMENT, NOT READY TO IMPLEMENT
 * -SIDD
 */



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

The minimum bracketing method is based on results from "Numerical Recipes 
in C, 2nd ed." and conforms to the authors' defintion of intellectual
property under their license, and as such does not conflict with
their copyright to their programs which execute similar algorithms.
*/

#include "nbody_priv.h"
#include "milkyway_util.h"
#include "milkyway_math.h"
#include "milkyway_lua.h"
#include "nbody_lua_types.h"
#include "nbody_isotropic.h"


/*Be Careful! this function returns the negative of the potential! this is the value of interest, psi*/
static inline real potential( real r, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
  real scaleRad1Cube = cube(scaleRad1); 
  real scaleRad2Cube = cube(scaleRad2);
  real potential= -1.0*(mass1/mw_sqrt(sqr(r) + sqr(scaleRad1)) +  Mass2/mw_sqrt(sqr(r) + sqr(scaleRad2)) );
  
  return (-1.0*potential);
}


static inline real density( real r, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
  real scaleRad1Cube = cube(scaleRad1); 
  real scaleRad2Cube = cube(scaleRad2);
  real density= (3/(4 M_pi)) (mass1/scaleRad1Cube *pow(1+ sqr(r)/sqr(scaleRad1), -2.5)  
			      + mass2/scaleRad2Cube *pow(1+ sqr(r)/sqr(scaleRad2), -2.5);
  
  return density;
}


static inline real fun(real ri, real mass1, real mass2, real scaleRad1, real scaleRad2, real upperlimit)
{
 real first_deriv_psi;
 real second_deriv_psi;
 real first_deriv_density;
 real second_deriv_density;
 real dsqden_dpsisq;/*second derivative of density with respect to -potential */
 real denominator; /*the demoninator of the distribution function: 1/sqrt(E-Psi)*/
 real func;
 real h;
 /*yes, this does in fact use a 5-point stencil*/
 first_deriv_psi=( potential(ri-2.0*h,mass1,mass2,scaleRad1,scaleRad2)-8.0*potential(ri-1.0*h,mass1,mass2,scaleRad1,scaleRad2)
		      -potential(ri+2.0*h,mass1,mass2,scaleRad1,scaleRad2)+8.0*potential(ri+1.0*h,mass1,mass2,scaleRad1,scaleRad2) ) /(12*h);
  
  first_deriv_density=( density(ri-2.0*h,mass1,mass2,scaleRad1,scaleRad2)-8.0*density(ri-1.0*h,mass1,mass2,scaleRad1,scaleRad2)
		      -density(ri+2.0*h,mass1,mass2,scaleRad1,scaleRad2)+8.0*density(ri+1.0*h,mass1,mass2,scaleRad1,scaleRad2) ) /(12*h);

/*yes, this also uses a five point stencil*/
  second_deriv_density=( -1.0*density(ri+2.0*h,mass1,mass2,scaleRad1,scaleRad2)+16.0*density(ri+1.0*h,mass1,mass2,scaleRad1,scaleRad2) -30.0*density(ri,mass1,mass2,scaleRad1,scaleRad2)
		      +16.0*density(ri-1.0*h,mass1,mass2,scaleRad1,scaleRad2)-1.0*density(ri-2.0*h,mass1,mass2,scaleRad1,scaleRad2) ) /(12*h*h);

  second_deriv_psi= ( -1.0*potential(ri+2.0*h,mass1,mass2,scaleRad1,scaleRad2)+16.0*potential(ri+1.0*h,mass1,mass2,scaleRad1,scaleRad2) -30.0*potential(ri,mass1,mass2,scaleRad1,scaleRad2)
		      +16.0*potential(ri-1.0*h,mass1,mass2,scaleRad1,scaleRad2)-1.0*potential(ri-2.0*h,mass1,mass2,scaleRad1,scaleRad2) ) /(12*h*h);

	/*
	 * Instead of calculating the second derivative of density with respect to -pot directly, 
	 * did product rule since both density and pot are functions of radius. 
	 */
  dsqden_dpsisq=second_deriv_density/ first_deriv_psi - first_deriv_density*second_deriv_psi/(mw_sqr(first_deriv_psi));
  demoninator= 1/mw_sqrt(upperlimit-potential(ri,mass1,mass2,scaleRad1,scaleRad2) );
  func= first_deriv_psi* dsqden_dpsisq *denominator;

  return func;
  
}
  
  
/*This is a guassian quadrature routine. It uses 1000 steps, so it should be quite accurate*/
static inline real gauss_quad( real lowerlimit, real upperlimit, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
  real Ng,hg,lowerg, upperg;
  real intv;
  real coef1,coef2;//parameters for gaussian quad
  real c1,c2,c3;
  real x1,x2,x3;
  real x1n,x2n,x3n;
  

  intv=0;//initial value of integral
  Ng=1001;
  hg=(b-a)/(Ng-1);

  lowerg=lowerlimit;
  upperg=lowerlimit+hg;
  coef2= (lowerg+upperg)/2;//initializes the first coeff to change the function limits
  coef1= (upperg-lowerg)/2;//initializes the second coeff to change the function limits
  c1=0.555555556;
  c2=0.888888889;
  c3=0.555555556;
  x1=-0.774596669;
  x2=0.000000000;
  x3=0.774596669;
  x1n=((coef1)*x1 +coef2);
  x2n=((coef1)*x2 +coef2);
  x3n=((coef1)*x3 +coef2);

  while (1)
  {

      //gauss quad
      intv= intv +(c1*fun(x1n, mass1, mass2, scaleRad1, scaleRad2, upperlimit)*coef1 +      
		    c2*fun(x2n, mass1, mass2, scaleRad1, scaleRad2, upperlimit)*coef1 + 
		    c3*fun(x3n, mass1, mass2, scaleRad1, scaleRad2, upperlimit)*coef1);

      lowerg=upperg;
      upperg= upperg+hg;
      coef2= (lowerg+ upperg)/2;//initializes the first coeff to change the function limits
      coef1= (upperg-lowerg)/2;
      x1n=((coef1)*x1 +coef2);
      x2n=((coef1)*x2 +coef2);
      x3n=((coef1)*x3 +coef2);

      if (lowerg>=upperlimit)//loop termination clause
        {break;}
  }

  return intv;
}

static inline real dist_fun(real mass1, real mass2, real scaleRad1, real scaleRad2, real lowerlimit, real upperlimit)
{
 real c= 1.0/(mw_sqrt(8)* cube(mw_pi));
 real distribution_function;
 distribution_function=c*gauss_quad(lowerlimit, upperlimit, mass1, mass2, scaleRad1, scaleRad2);
  
  return distribution_function;
}

/* pickshell: pick a random point on a sphere of specified radius. 
*
*   Changed this section to be in compliance with NEMO's initialization technique.
*   Instead of assigning it a point on sphere, assigns angles. Allows for non-circular orbits.
*
*/
/*NEED TO PERHAPS CHANGE*/
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

/*NEED TO CHANGE*/
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


/*NEED TO CHANGE
 *- changed the first profile
 * need to specify p_crit
 */
static inline real profile(real r, real mass1, real mass2, real scale1, real scale2, real p_0)
  {  
    //real p_crit= 3*H*H/(8*M_pi*G);//have to input values
    real prof = (  r*r *( p_0*mass1*scale1/(r)* mw_pow(1.0 + r/ scale1,-2.0) 
				 + (mass2/mw_pow(scale2,3.0)) * mw_pow(1 + mw_pow(r,2.0) / mw_pow(scale2,2.0),-2.5) )   );
    return (real) (-prof);
  }

  


/*NEED TO CHANGE*/
real computeRhoMax(real mass1, real mass2, real scale1, real scale2)
{
  real result =  inverseParabolicInterpolateIsotropic(0.61803399, 0, 5.0 * (scale1 + scale2), 
				       scale2, scale1, scale2, mass1, mass2, 1e-4);  
  mw_printf("RHO MAX IS %10.5f\n",result);
  return result;
}

/*NEED TO CHANGE
 * changed thefirst profile. 
 * need to specify p_crit and deltac
 */
static inline real isotropicRandomR(dsfmt_t* dsfmtState, real scaleRad1, real scaleRad2,
				    real Mass1, real Mass2,real max)
{

  real scaleRad2Cube = cube(scaleRad2);
  //real p_crit= 3*H*H/(8*M_pi*6.67e-11);//have to input values
  real delta=1.0;//have to input value
  mwbool GOOD_RADIUS = 0;

  real r;
  real u, val;

  while (GOOD_RADIUS != 1)
    {
      r = (real)mwXrandom(dsfmtState,0.0, 5.0 * (scaleRad1 + scaleRad2));
      u = (real)mwXrandom(dsfmtState,0.0,1.0);

      val = r*r * (  Mass1*scaleRad1/(r)* mw_pow(1.0 + r/ scaleRad1,-2.0) +
					( 3.0/(4.0 *M_PI)*Mass2/scaleRad2Cube * mw_pow(1.0 + sqr(r)/sqr(scaleRad2),-2.5) )  );

      if (val/max > u)
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

/*NEED TO CHANGE*/
static inline real isotropicRandomV(dsfmt_t* dsfmtState,real r, real scaleRad1, real scaleRad2, real Mass1, real Mass2)
{

  real GMsolar = 1.327e20; //SI UNITS m^3 / sec^2
  scaleRad1 *= 3.086e19; //meters
  scaleRad2 *= 3.086e19;  
  r *= 3.086e19;
  real val;
  real x;
  x= plummerSelectFromG(dsfmtState);
/*this calculates it in m/s. return val is converted to km/s thus mult by 0.001*/
  val = x * M_SQRT2* mw_sqrt(GMsolar *Mass1/mw_sqrt(sqr(r) + sqr(scaleRad1))
			       + GMsolar * Mass2/mw_sqrt(sqr(r) + sqr(scaleRad2)));

  return 0.001* val; //km/s
}

/*NEED TO CHANGE*/
static inline mwvector isotropicBodyPosition(dsfmt_t* dsfmtState, mwvector rshift,  real r)
{
    mwvector pos;

    pos = pickShell(dsfmtState,  r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */

    return pos;
}

/*NEED TO CHANGE*/
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

    /*NEED TO CHANGE*/velScale = 1000;// Conversion from km/s

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);
    /*NEED TO CHANGE*/ real RHO_MAX = computeRhoMax((real)mass1, (real)mass2, (real)radiusScale1, (real)radiusScale2);
    mw_printf("%10.5f",RHO_MAX);

    for (i = 0; i < nbody; ++i)
    {
        do
        {
         /*NEED TO CHANGE*/ r = isotropicRandomR(prng, radiusScale1, radiusScale2, mass1, mass2, RHO_MAX);
	           /* FIXME: We should avoid the divide by 0.0 by multiplying
             * the original random number by 0.9999.. but I'm too lazy
             * to change the tests. Same with other models */
        }
        while (isinf(r));

        /*NEED TO CHANGE*/b.bodynode.pos = isotropicBodyPosition(prng, rShift, r);

        /*NEED TO CHANGE*/b.vel = isotropicBodyVelocity(prng, r, vShift, velScale, radiusScale1, radiusScale2, mass1, mass2);

        assert(nbPositionValid(b.bodynode.pos));

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

/*NEED TO CHANGE LUA FILE*/
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

/*NEED TO CHANGE*/
void registerGenerateIsotropic(lua_State* luaSt)
{
    /*NEED TO CHANGE*/lua_register(luaSt, "generateIsotropic", nbGenerateIsotropic);
}


