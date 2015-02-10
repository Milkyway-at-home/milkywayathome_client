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
#include <time.h>
#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */




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
  real density_result= (3.0/(4.0*(M_PI)))*( (mass1/scaleRad1Cube) *mw_pow(1.0 + sqr(r)/sqr(scaleRad1), -2.5)
						  + (mass2/scaleRad2Cube) *mw_pow(1.0 + sqr(r)/sqr(scaleRad2), -2.5));
  
  return density_result;
}

/*BE CAREFUL! this function returns the mass enclosed in a single plummer sphere!*/
static inline real mass_en( real r, real mass, real scaleRad)
{
  real mass_enclosed= mass* cube(r)*pow( (sqr(r)+ sqr(scaleRad)), -1.5);
  
  return mass_enclosed;
}


static inline real fun(real ri, real mass1, real mass2, real scaleRad1, real scaleRad2, real energy)
{
 real first_deriv_psi;
 real second_deriv_psi;
 real first_deriv_density;
 real second_deriv_density;
 real dsqden_dpsisq;/*second derivative of density with respect to -potential (psi) */
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
  denominator= 1.0/mw_sqrt(mw_fabs(energy-potential(ri,mass1,mass2,scaleRad1,scaleRad2) ));
  func= first_deriv_psi* dsqden_dpsisq *denominator;

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
      

      if (lowerg>=energy)//loop termination clause
        {break;}
  }
  return intv;
}

 /*This returns the value of the distribution function for a given energy*/
static inline real dist_fun(real r, real mass1, real mass2, real scaleRad1, real scaleRad2, real energy)
{
 real c= 1.0/(mw_sqrt(8)* sqr(M_PI));
 real distribution_function;
/*This calls guassian quad to integrate the function for a given energy*/
 distribution_function=c*gauss_quad(energy, mass1, mass2, scaleRad1, scaleRad2);
  
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

    /*legacy code*/
//     rsq = mw_sqrv(vec);             /* compute radius squared */
//     rsc = rad / mw_sqrt(rsq);       /* compute scaling factor */
//     mw_incmulvs(vec, rsc);          /* rescale to radius given */

    return vec;
}

 

/*this serves the max finding routine*/
static inline real profile2(real r, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
    real result =  r*r*density( r, mass1, mass2, scaleRad1, scaleRad2);  
  return result;
}
  
  static inline real profile(real v, real mass1, real mass2, real r, real scaleRad1, real scaleRad2 )
{
  real energy= potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*v*v;
  real result =  v*v*dist_fun( r, mass1, mass2, scaleRad1, scaleRad2, energy);  
  return -result;
}

/*this is a maxfinding routine to find the maximum of the density.
 * It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
 */

static inline real rhomax_finder( real a, real b, real c, real scaleRad1, real scaleRad2, real mass1,real mass2)
{
  real tolerance= 1e-4;
  real RATIO = 0.61803399;
  real RATIO_COMPLEMENT = 1 - RATIO;
  int counter=0;
  
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
      counter++;
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
	if(counter>20){break;}
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

static inline real distmax_finder( real a, real b, real c, real r, real scaleRad1, real scaleRad2, real mass1,real mass2)
{
  real tolerance= 1e-2;
  real RATIO = 0.61803399;
  real RATIO_COMPLEMENT = 1 - RATIO;
  int counter=0;
  
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

  profile_x1 = (real)profile(x1,mass1,mass2,r, scaleRad1,scaleRad2);
  profile_x2 = (real)profile(x2,mass1,mass2,r, scaleRad1,scaleRad2);
  
  while (mw_fabs(x3 - x0) > (tolerance * (mw_fabs(x1) + mw_fabs(x2)) ) )
    {
      counter++;
      if (profile_x2 < profile_x1)
	{
	  x0 = x1;
	  x1 = x2;
	  x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
	  profile_x1 = (real)profile_x2;
	  profile_x2 = (real)profile(x2, mass1,mass2,r, scaleRad1,scaleRad2);
	}
      else
	{
	  x3 = x2;
	  x2 = x1;
	  x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
	  profile_x2 = (real)profile_x1;
	  profile_x1 = (real)profile(x1,mass1,mass2,r, scaleRad1,scaleRad2);
	}
	
      if(counter>10){break;}
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
//       mw_printf(" \n  %f      %f ",r, u);
      val = r*r * density(r,  mass1,  mass2,  scaleRad1,  scaleRad2);
  
      if (val/rho_max > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  return r;
}



static inline real vel_mag(dsfmt_t* dsfmtState,real r, real mass1, real mass2, real scaleRad1, real scaleRad2,  real part_mass)
{
  
  /*
   * WE TOOK IN 
   * MASS IN SIMULATION UNITS, WHICH HAVE THE UNITS OF KPC^3/GY^2 
   * LENGTH IN KPC
   * AND TIME IN GY
   * THEREFORE, THE velocities ARE OUTPUTING IN KPC/GY
   * THIS IS EQUAL TO 0.977813107 KM/S
   * 
   */
  
  
      
//   real GMsolar =222288.47; //convert from simulation to solar masses
//   scaleRad1 *= 1000; //pc
//   scaleRad2 *= 1000;  
//   r *= 1000;
//   mass1 *=GMsolar;
//   mass2 *=GMsolar;
//   part_mass*=GMsolar;
  
  real val,v,u,d;
  real energy;

  //real v_esc= mw_sqrt( mw_fabs(2.0*potential(r, mass1, mass2, scaleRad1,scaleRad2))); 
  real v_esc= mw_sqrt( mw_fabs(2.0* (mass1+mass2)/r));
  real dist_max=distmax_finder( 0.0, .5*v_esc, v_esc, r, scaleRad1,  scaleRad2, mass1, mass2);
    while (1)
    {
      v = (real)mwXrandom(dsfmtState,0.0, v_esc);
      u = (real)mwXrandom(dsfmtState,0.0,1.0);
      
      energy= potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*v*v;
      
      d=dist_fun( r,  mass1,  mass2,  scaleRad1,  scaleRad2, energy);
      val =v*v* d;
      if (mw_fabs( val/dist_max) > u)
      {
       	break;
      }
    }
  
  
  v*=0.977813107;//changing from kpc/gy to km/s
  return v; //km/s
}


static inline mwvector r_vec(dsfmt_t* dsfmtState, mwvector rshift,  real r)
{
    mwvector pos;

    pos = angles(dsfmtState,  r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */

    return pos;
}


static inline mwvector vel_vec(dsfmt_t* dsfmtState, mwvector vshift,real v)
{
    mwvector vel;

    vel = angles(dsfmtState, v);   /* pick scaled velocity */
    mw_incaddv(vel, vshift);                /* move the velocity */

    return vel;
}

// static inline void dist_func_plot(real r, real mass1, real mass2, real scaleRad1, real scaleRad2,real part_mass)
// {
//   
//   real energy,d;
//   real v=0.0;
//   real v_esc= mw_sqrt( mw_fabs(2.0* (mass1+mass2)/r));  
//   FILE *fp;
//   fp= fopen("dist.dat", "w");  
//   while(1)
//   {
//     energy= part_mass*potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*part_mass*v*v;
//     d=dist_fun(r,mass1,mass2, scaleRad1, scaleRad2, energy);
//     
//     fprintf(fp, "%.16f \t %.16f %.16f %.16f \n", v, v*v*d,r, energy);
//     v+=0.001;
//     if(v>=v_esc){break;}
//   }
//   fclose(fp);
// }

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
    real r, v;
    real mass_en1, mass_en2; //mass enclosed within predetermined r
    real mass = mass1 + mass2;
    
    memset(&b, 0, sizeof(b));
    
    real rho_max=-rhomax_finder(0,radiusScale2, (radiusScale1 + radiusScale2), radiusScale1, radiusScale2, mass1, mass2);

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */
    real p_mass=mass / nbody;
    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);	
    real all_rs[nbody];
    real all_vs[nbody];


     
      for (i = 0; i < nbody; i++)
      {
// 	 mw_printf(" \r initalizing particle %i. ",i);
	  do
	  {
	    r= r_mag(prng, mass1, mass2, radiusScale1, radiusScale2, rho_max);
	    /*to ensure that r is finite and nonzero*/
	    if(isinf(r)==FALSE && r!=0.0){break;}
	  }
	  while (1);
	  all_rs[i]=r;

	  /*this calculates the mass enclosed in each sphere. 
	  * velocity is determined by mass enclosed at that r not by the total mass of the system. 
	  */
	  mass_en1= mass_en(r, mass1, radiusScale1);
	  mass_en2= mass_en(r, mass2, radiusScale2);
	  
	  all_vs[i] = vel_mag(prng, r, mass_en1, mass_en2, radiusScale1, radiusScale2, p_mass);
      }
    
    for(i=0;i<nbody;i++)
    {
      r=all_rs[i];
      v=all_vs[i];
      b.bodynode.pos = r_vec(prng, rShift, r);
      b.vel = vel_vec(prng,  vShift,v);
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


