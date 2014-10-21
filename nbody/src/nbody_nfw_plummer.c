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
  real potential= -1.0*(mass1/mw_sqrt(sqr(r) + sqr(scaleRad1)) +  mass2/mw_sqrt(sqr(r) + sqr(scaleRad2)) );
  
  return (-1.0*potential);
}

/*this is the density distribution function. Returns the density at a given radius.*/
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
 real h=0.01; /*This value is not completely arbitrary. Generally, lower the value of h the better. 
 For the five point stencil, values lower than .001 ran into roundoff error.
 0.01 is a safe bet, with a calculation error of order 10^-10.*/
 
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
static inline real gauss_quad(  real upperlimit, real mass1, real mass2, real scaleRad1, real scaleRad2)
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
/*I have set the lower limit to be zero. '
 * This is in the definition of the distribution function. 
 * If this is used for integrating other things, this will need to be changed.*/
  lowerg=0.0;
  upperg=lowerg+hg;
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

 /*This returns the value of the distribution function for a given energy*/
static inline real dist_fun(real r, real mass1, real mass2, real scaleRad1, real scaleRad2, real upperlimit)
{
 real c= 1.0/(mw_sqrt(8)* sqr(mw_pi));
 real distribution_function;
/*This calls guassian quad to integrate the function for a given energy*/
 distribution_function=c*gauss_quad(upperlimit, mass1, mass2, scaleRad1, scaleRad2);
  
  return distribution_function;
}

/* assigns angles. Allows for non-circular orbits.*/
static inline mwvector angles(dsfmt_t* dsfmtState, real rad)
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

/*NEED TO CHANGE*/

static inline real profile(real r, real mass1, real mass2, real scale1, real scale2, real p_0)
  {  
    real prof = (  r*r *density( r,  mass1,  mass2,  scaleRad1,  scaleRad2));
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


static inline real r_mag(dsfmt_t* dsfmtState, real mass1, real mass2, real scaleRad1, real scaleRad2,real rho_max)
{

  mwbool GOOD_RADIUS = 0;

  real r;
  real u, val;

  while (GOOD_RADIUS != 1)
    {
      r = (real)mwXrandom(dsfmtState,0.0, 5.0 * (scaleRad1 + scaleRad2));
      u = (real)mwXrandom(dsfmtState,0.0,1.0);

      val = r*r * density(r,  mass1,  mass2,  scaleRad1,  scaleRad2);

      if (val/rho_max > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  return r;
}


/*NEED TO CHANGE*/
static inline real vel_mag(dsfmt_t* dsfmtState,real r, real mass1, real mass2, real scaleRad1, real scaleRad2,  real dist_max)
{

  real GMsolar = 1.327e20; //SI UNITS m^3 / sec^2
  scaleRad1 *= 3.086e19; //meters
  scaleRad2 *= 3.086e19;  
  r *= 3.086e19;
  real val,v;
  /*there is a -1 there because potential actually returns the neg of potential*/
  real v_esc= mw_sqrt( -2.0*potential(r, mass1, mass2, radiusScale1,radiusScale2)); 
  real dist_max = max_finder(r, mass1, mass2, radiusScale1,radiusScale2);
/*this calculates it in m/s. return val is converted to km/s thus mult by 0.001*/


    while (GOOD_RADIUS != 1)
    {
      v = (real)mwXrandom(dsfmtState,0.0, v_esc);
      u = (real)mwXrandom(dsfmtState,0.0,1.0);

      val = 4*mw_pi*v*v* dist_fun( r,  mass1,  mass2,  scaleRad1,  scaleRad2);

      if (val/dist_max > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  
  
  return 0.001* val; //km/s
}


static inline mwvector r_vec(dsfmt_t* dsfmtState, mwvector rshift,  real r)
{
    mwvector pos;

    pos = angles(dsfmtState,  r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */

    return pos;
}


static inline mwvector vel_vec(dsfmt_t* dsfmtState, mwvector vshift, real vsc,real r, real mass1, real mass2,
			       real scaleRad1, real scaleRad2)
{
    mwvector vel;
    real v;

    v = vel_mag(dsfmtState, r, mass1, mass2, scaleRad1, scaleRad2);
    vel = angles(dsfmtState, vsc * v);   /* pick scaled velocity */
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
    real RHO_MAX = computeRhoMax((real)mass1, (real)mass2, (real)radiusScale1, (real)radiusScale2);
    mw_printf("%10.5f",RHO_MAX);

    for (i = 0; i < nbody; ++i)
    {
        do
        {
         r = r_mag(prng, mass1, mass2, radiusScale1, radiusScale2, RHO_MAX);
	 
	          /* FIXME: We should avoid the divide by 0.0 by multiplying
             * the original random number by 0.9999.. but I'm too lazy
             * to change the tests. Same with other models */
        }
        while (isinf(r));

        /*NEED TO CHANGE*/b.bodynode.pos = r_vec(prng, rShift, r);

        /*NEED TO CHANGE*/b.vel = vel_vec(prng,  vShift, velScale,r, mass1, mass2, radiusScale1, radiusScale2);

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

