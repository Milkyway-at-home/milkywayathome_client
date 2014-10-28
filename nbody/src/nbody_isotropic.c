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

The minimum bracketing method is based on results from "Numerical Recipes,
3rd ed." and conforms to the authors' defintion of intellectual
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
  real potential_result= -1.0*(mass1/mw_sqrt(sqr(r) + sqr(scaleRad1)) +  mass2/mw_sqrt(sqr(r) + sqr(scaleRad2)) );
  
  return (-1.0*potential_result);
}

/*this is the density distribution function. Returns the density at a given radius.*/
static inline real density( real r, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
  real scaleRad1Cube = cube(scaleRad1); 
  real scaleRad2Cube = cube(scaleRad2);
  real density_result= (3/(4*(M_PI)))*(mass1/scaleRad1Cube *pow(1+ sqr(r)/sqr(scaleRad1), -2.5)+ mass2/scaleRad2Cube *pow(1+ sqr(r)/sqr(scaleRad2), -2.5));
  
  return density_result;
}


static inline real fun(real ri, real mass1, real mass2, real scaleRad1, real scaleRad2, real energy)
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
  dsqden_dpsisq=second_deriv_density/ first_deriv_psi - first_deriv_density*second_deriv_psi/(sqr(first_deriv_psi));
  denominator= 1/mw_sqrt(mw_fabs(-potential(ri,mass1,mass2,scaleRad1,scaleRad2) ));
  func= first_deriv_psi* dsqden_dpsisq *denominator;
  //mw_printf("function value= %f \n",func);
  return func;
  
}
  
  
/*This is a guassian quadrature routine. It uses 1000 steps, so it should be quite accurate*/
static inline real gauss_quad(  real energy, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
  real Ng,hg,lowerg, upperg;
  real intv;
  real coef1,coef2;//parameters for gaussian quad
  real c1,c2,c3;
  real x1,x2,x3;
  real x1n,x2n,x3n;
  real a=0.0;
  real b=energy;

  intv=0;//initial value of integral
  Ng=1001;
  hg=(b-a)/(Ng-1);
/*I have set the lower limit to be zero. '
 * This is in the definition of the distribution function. 
 * If this is used for integrating other things, this will need to be changed.*/
  lowerg=0.0;
  upperg=lowerg+hg;
  
//   mw_printf("upper= %f  energy= %f \n", upperg, energy);
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

//   mw_printf("       integrating...\n");
  while (1)
  {
      
      //gauss quad
      intv= intv +(c1*fun(x1n, mass1, mass2, scaleRad1, scaleRad2, energy)*coef1 +      
		    c2*fun(x2n, mass1, mass2, scaleRad1, scaleRad2, energy)*coef1 + 
		    c3*fun(x3n, mass1, mass2, scaleRad1, scaleRad2, energy)*coef1);

      lowerg=upperg;
      upperg= upperg+hg;
      coef2= (lowerg+ upperg)/2;//initializes the first coeff to change the function limits
      coef1= (upperg-lowerg)/2;
      x1n=((coef1)*x1 +coef2);
      x2n=((coef1)*x2 +coef2);
      x3n=((coef1)*x3 +coef2);
      
//    mw_printf("lower= %f  energy= %f", lowerg, energy);

      if (lowerg>=energy)//loop termination clause
        {break;}
  }
//       mw_printf("       done integrating \n");
  return intv;
}

 /*This returns the value of the distribution function for a given energy*/
static inline real dist_fun(real r, real mass1, real mass2, real scaleRad1, real scaleRad2, real energy)
{
 real c= 1.0/(mw_sqrt(8)* sqr(M_PI));
 real distribution_function;
/*This calls guassian quad to integrate the function for a given energy*/
//   mw_printf("      entering integrating routine... \n");
 distribution_function=c*gauss_quad(energy, mass1, mass2, scaleRad1, scaleRad2);
//  mw_printf("      finished integrating. \n");
  
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

/*this serves the max finding routine*/
static inline real profile(real v, real r, real mass1, real mass2, real scaleRad1, real scaleRad2, real part_mass)
{
    real energy= potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*part_mass*v*v;
//     mw_printf(" v= %f \n",  v);
//     mw_printf("     getting distribution function for initial profiles...\n");
    real result =  dist_fun( r, mass1, mass2, scaleRad1, scaleRad2, energy);  
//     mw_printf("     retreived distribution function for initial profiles \n");
  return result;
}


  
/*this is a maxfinding routine to find the maximum of the distribution function for a 
 * given radius. It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
 */
static inline real distmax_finder( real a, real b, real c,
		     real r, real scaleRad1, real scaleRad2, real mass1,real mass2, real part_mass)
{
//   mw_printf("   performing max finding routine: a= %f, b= %f , c= %f \n",a, b, c);
  int counter=0;
  int interval=1000;
  real tolerance= 1e-4;
  real RATIO = 0.61803399;
  real RATIO_COMPLEMENT = 1 - RATIO;
  
  real profile_x1,profile_x2,x0,x1,x2,x3;
  x0 = a;
  x3 = c;
  
  if (mw_fabs(b - c) > mw_fabs(b - a))
    {
      x1 = b;
      x2 = b + (RATIO_COMPLEMENT * (c - b)); 
    }
  else
    {
      x2 = b;
      x1 = b - (RATIO_COMPLEMENT * (b - a));
    }
//   mw_printf("    getting initial profiles...\n");
  profile_x1 = (real)profile(x1,r,mass1,mass2,scaleRad1,scaleRad2, part_mass);
  profile_x2 = (real)profile(x2,r,mass1,mass2,scaleRad1,scaleRad2, part_mass);
//   mw_printf("    profiles retreived \n");
  
  while (mw_fabs(x3 - x0) > (tolerance * (mw_fabs(x1) + mw_fabs(x2)) ) )
    {
//       if ((counter%interval)==0){mw_printf("counter = %i \n",counter);}
//      counter++;
      if (profile_x2 < profile_x1)
	{
	  x0 = x1;
	  x1 = x2;
	  x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
	  profile_x1 = (real)profile_x2;
	  profile_x2 = (real)profile(x2, r, mass1,mass2,scaleRad1,scaleRad2, part_mass);
	}
      else
	{
	  x3 = x2;
	  x2 = x1;
	  x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
	  profile_x2 = (real)profile_x1;
	  profile_x1 = (real)profile(x1,r,mass1,mass2,scaleRad1,scaleRad2, part_mass);
	}
    }
//     mw_printf("counter = %i \n",counter);
//  mw_printf("   finished max finding routine");
  if (profile_x1 < profile_x2)
    {
      return (-profile_x1);
    }
  else
    {
      return (-profile_x2);
    }
}

/*this serves the max finding routine*/
static inline real profile2(real r, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
    real result =  r*r*density( r, mass1, mass2, scaleRad1, scaleRad2);  
  return result;
}
  
/*this is a maxfinding routine to find the maximum of the density.
 * It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
 */
static inline real rhomax_finder( real a, real b, real c, real scaleRad1, real scaleRad2, real mass1,real mass2)
{
  real tolerance= 1e-4;
  real RATIO = 0.61803399;
  real RATIO_COMPLEMENT = 1 - RATIO;
//   int counter=0;
  
  real profile_x1,profile_x2,x0,x1,x2,x3;
  x0 = a;
  x3 = c;
  
  if (mw_fabs(b - c) > mw_fabs(b - a))
    {
      x1 = b;
      x2 = b + (RATIO_COMPLEMENT * (c - b)); 
    }
  else
    {
      x2 = b;
      x1 = b - (RATIO_COMPLEMENT * (b - a));
    }

  profile_x1 = (real)profile2(x1,mass1,mass2,scaleRad1,scaleRad2);
  profile_x2 = (real)profile2(x2,mass1,mass2,scaleRad1,scaleRad2);
  
  while (mw_fabs(x3 - x0) > (tolerance * (mw_fabs(x1) + mw_fabs(x2)) ) )
    {
//       counter++;
      if (profile_x2 < profile_x1)
	{
	  x0 = x1;
	  x1 = x2;
	  x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
	  profile_x1 = (real)profile_x2;
	  profile_x2 = (real)profile2(x2, mass1,mass2,scaleRad1,scaleRad2);
	}
      else
	{
	  x3 = x2;
	  x2 = x1;
	  x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
	  profile_x2 = (real)profile_x1;
	  profile_x1 = (real)profile2(x1,mass1,mass2,scaleRad1,scaleRad2);
	}
    }
//   mw_printf("counter = %i \n",counter);
  if (profile_x1 < profile_x2)
    {
      return (-profile_x1);
    }
  else
    {
      return (-profile_x2);
    }
}

static inline real r_mag(dsfmt_t* dsfmtState, real mass1, real mass2, real scaleRad1, real scaleRad2, real rho_max)
{

  mwbool GOOD_RADIUS = 0;

  real r;
  real u, val;

  while (GOOD_RADIUS != 1)
    {
      r = (real)mwXrandom(dsfmtState,0.0, 5.0 * (scaleRad1 + scaleRad2));
      u = (real)mwXrandom(dsfmtState,0.0,1.0);

//       mw_printf("getting density at r...");
      val = r*r * density(r,  mass1,  mass2,  scaleRad1,  scaleRad2);
//       mw_printf("done. val= %f \n",val/rho_max);
  
      if (val/rho_max > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  return r;
}


/*NEED TO CHANGE*/
static inline real vel_mag(dsfmt_t* dsfmtState,real r, real mass1, real mass2, real scaleRad1, real scaleRad2,  real part_mass)
{
  mwbool GOOD_VAL= 0;
//   real GMsolar = 1.327e20; //SI UNITS m^3 / sec^2
//   scaleRad1 *= 3.086e19; //meters
//   scaleRad2 *= 3.086e19;  
//   r *= 3.086e19;
  real val,v,u;
  real energy;
  /*there is a -1 there because potential actually returns the neg of potential*/
  real v_esc= mw_sqrt( mw_fabs(2.0*potential(r, mass1, mass2, scaleRad1,scaleRad2))); 
//   mw_printf("   vesc= %f \n", v_esc);
//   mw_printf("   getting distribution_function max...\n");
  /*This is:	   distmax_finder(a,   b,           c,                       r, ... particle mass*/
  real dist_max = distmax_finder(0.0, 0.5*v_esc, 1.0 * v_esc, r, scaleRad1, scaleRad2, mass1, mass2, part_mass);
//   mw_printf("   dist_max=%f \n", dist_max);
  
/*this calculates it in m/s. return val is converted to km/s thus mult by 0.001*/


    while (GOOD_VAL != 1)
    {
      v = (real)mwXrandom(dsfmtState,0.0, v_esc);
      u = (real)mwXrandom(dsfmtState,0.0,1.0);
      
//       mw_printf("   getting energy for distribution_function... \n");
      energy= potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*part_mass*v*v;
//       mw_printf("   done. energy= %f \n",energy);
      
//       mw_printf("   getting value...");
      val = 4*M_PI*v*v* dist_fun( r,  mass1,  mass2,  scaleRad1,  scaleRad2, energy);
//       mw_printf("   done. val= %f , val/dist_max= %f, u= %f \n" ,val, val/dist_max, u);
      if (mw_fabs( val/dist_max) > u)
      {
// 	mw_printf("test complete...val/dist_max= %f, u= %f \n" , val/dist_max, u);
       	GOOD_VAL = 1;
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
			       real scaleRad1, real scaleRad2, real part_mass)
{
    mwvector vel;
    real v;
    
//     mw_printf("  getting velocity magnitude...\n");
    v = vel_mag(dsfmtState, r, mass1, mass2, scaleRad1, scaleRad2, part_mass);
//     mw_printf("  done. vel= %f \n", vel_mag);
    
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
    
//     mw_printf("getting rhomax...");
    real rho_max=-rhomax_finder(0,radiusScale2, 5.0 * (radiusScale1 + radiusScale2), radiusScale1, radiusScale2, mass1, mass2);
//     mw_printf("done. rhomax= %f \n", rho_max);
    
    /*NEED TO CHANGE*/velScale = 1000;// Conversion from km/s

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < nbody; ++i)
    {
//       mw_printf("initalizing particle %i \n",i);
        do
        {
// 	  mw_printf(" getting radius for particle %i...", i);
         r = r_mag(prng, mass1, mass2, radiusScale1, radiusScale2, rho_max);
	 
	          /* FIXME: We should avoid the divide by 0.0 by multiplying
             * the original random number by 0.9999.. but I'm too lazy
             * to change the tests. Same with other models */
        }
        while (isinf(r));
// 	  mw_printf(" done, r= %f \n",r);
// 	mw_printf(" getting radius vector...");
        /*NEED TO CHANGE*/b.bodynode.pos = r_vec(prng, rShift, r);
// 	mw_printf(" done \n");
	
// 	mw_printf(" getting velocity vector...\n");
        /*NEED TO CHANGE*/b.vel = vel_vec(prng,  vShift, velScale,r, mass1, mass2, radiusScale1, radiusScale2, b.bodynode.mass);
// 	mw_printf(" done \n");
        assert(nbPositionValid(b.bodynode.pos));

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
	
// 	mw_printf("finished particle %i.", i);
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


