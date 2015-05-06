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

/*Note: minusfivehalves(x) raises to x^-5/2 power and minushalf(x) is x^-1/2*/

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
  /*this weird sqrt(fifth(x) notation is used because it was determined that pow() ate up a lot of comp time*/
  real density_result= (3.0/(4.0*(M_PI)))*( (mass1/scaleRad1Cube) * (minusfivehalves( (1.0 + sqr(r)/sqr(scaleRad1)) )  )
						  + (mass2/scaleRad2Cube) *(minusfivehalves( (1.0 + sqr(r)/sqr(scaleRad2))  ) ) );
  
  return density_result;
}


/*BE CAREFUL! this function returns the mass enclosed in a single plummer sphere!*/
static inline real mass_en( real r, real mass, real scaleRad)
{
  real mass_enclosed= mass* cube(r)* minusthreehalves( (sqr(r)+ sqr(scaleRad) ) ) ;
  
  return mass_enclosed;
}



static  real fun(real ri, real mass1, real mass2, real scaleRad1, real scaleRad2, real energy)
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
  denominator= minushalf( mw_fabs(energy-potential(ri,mass1,mass2,scaleRad1,scaleRad2) ));
  func= first_deriv_psi* dsqden_dpsisq *denominator;

  return func;
  
}


/*This is a guassian quadrature routine. It uses 1000 steps, so it should be quite accurate*/
static real gauss_quad(real upper, real energy, real mass1, real mass2, real scaleRad1, real scaleRad2)
{
  real Ng,hg,lowerg, upperg;
  real intv;
  real coef1,coef2;//parameters for gaussian quad
  real c1,c2,c3;
  real x1,x2,x3;
  real x1n,x2n,x3n;
  
  //this should be from infinity. But the dis func should be negligble here.
  real a=10.0*(scaleRad1+scaleRad2);
  real b=upper;
  

  intv=0;//initial value of integral
  Ng=1001.0;//integral resolution
  hg=fabs(b-a)/(Ng-1.0);
/*I have set the lower limit to be zero. '
 * This is in the definition of the distribution function. 
 * If this is used for integrating other things, this will need to be changed.*/
  lowerg=0.0;
  upperg=lowerg+hg;
  

  coef2= (lowerg+upperg)/2.0;//initializes the first coeff to change the function limits
  coef1= (upperg-lowerg)/2.0;//initializes the second coeff to change the function limits
  c1=0.555555556;
  c2=0.888888889;
  c3=0.555555556;
  x1=-0.774596669;
  x2=0.000000000;
  x3=0.774596669;
  x1n=((coef1)*x1 +coef2);
  x2n=((coef1)*x2 +coef2);
  x3n=((coef1)*x3 +coef2);
int counter=0;
  while (1)
  {

      //gauss quad
      intv= intv +(c1*fun(x1n, mass1, mass2, scaleRad1, scaleRad2, energy)*coef1 +      
		    c2*fun(x2n, mass1, mass2, scaleRad1, scaleRad2, energy)*coef1 + 
		    c3*fun(x3n, mass1, mass2, scaleRad1, scaleRad2, energy)*coef1);

      lowerg=upperg;
      upperg= upperg+hg;
      coef2= (lowerg+ upperg)/2.0;//initializes the first coeff to change the function limits
      coef1= (upperg-lowerg)/2.0;
      
      x1n=((coef1)*x1 +coef2);
      x2n=((coef1)*x2 +coef2);
      x3n=((coef1)*x3 +coef2);

      
      if (lowerg>=b)//loop termination clause
        {break;}
        counter++;
  }
  return intv;
}

/* A root finding function used to invert probability functions allows for direct calculation of positions instead of sampling */
static real findRoot(real (*rootFunc)(real, real*), real* rootFuncParams, real funcValue, real lowBound, real upperBound, dsfmt_t* dsfmtState)
{
  if(rootFuncParams == NULL || rootFunc == NULL)
  {exit(-1);}
  
  real xo;
  real fun, m, b, temp, tempf;
  real forward, backward;
  real h=0.001;
  int N=5;
  real roots[5]={0,0,0,0,0};
  int counter=0.0;
  int counter2=0.0;
  int found_root=0;
  for(int i=0; i<N; i++)
  {
    xo= mwXrandom(dsfmtState,0.0,1.0)*(upperBound-lowBound) +lowBound;
    counter=0;
    while(1)
    {
      counter2=0.0;
      fun=(*rootFunc)(xo, rootFuncParams)- funcValue;

      if(fabs(fun)<1e-2)
      {
	found_root++;
	break;
      }
      if(counter>=1000)
      {
	break;
      }
       counter++;
       
      forward=(*rootFunc)(xo+h, rootFuncParams)- funcValue;
      backward=(*rootFunc)(xo-h, rootFuncParams)- funcValue;
      m= (forward-backward)/(2.0*h);
      b= fun-(m*xo);
      
//       backtracking:
      temp=-1.0*b/m;
      tempf=(*rootFunc)(temp, rootFuncParams)-funcValue;
      while(fabs(tempf)>fabs(fun))
      {
	counter2++;
	if(counter2>=1000)
	{
	  break;
	}
	temp=temp/2.0;
	tempf=(*rootFunc)(temp, rootFuncParams)-funcValue;
      }
      xo=temp;
    }
    roots[i]=xo;
  }
  
  if(found_root!=0)
  {
    real root=roots[(int)(mwXrandom(dsfmtState,0.0,.999999)*(int)5)];
    return root;
  }
  else
  {return 0.0;}
  
}


 /*This returns the value of the distribution function for a given energy*/
static real dist_fun(real mass1, real mass2, real scaleRad1, real scaleRad2, real energy, real upper)
{
 real c= inv( (mw_sqrt(8)* sqr(M_PI)) );
 real distribution_function;
 
/*This calls guassian quad to integrate the function for a given energy*/
 distribution_function=c*gauss_quad(upper,energy, mass1, mass2, scaleRad1, scaleRad2);
  
  return distribution_function;
}

static inline real profile_deriv_rho(real r, real * args)
{
  real mass1 = args[0];
  real mass2 = args[1];
  real scaleRad1 = args[2];
  real scaleRad2 = args[3];
  real h=0.001;
  real forward=sqr(r+h)*density( r+h, mass1, mass2, scaleRad1, scaleRad2);
  real backward=sqr(r-h)*density( r-h, mass1, mass2, scaleRad1, scaleRad2); 
  real deriv= (forward-backward)/(2.0*h);
  return deriv;
}

static inline real profile_psi(real r, real * args)
{
  real mass1 = args[0];
  real mass2 = args[1];
  real scaleRad1 = args[2];
  real scaleRad2 = args[3];
  real psi=potential(r, mass1, mass2,scaleRad1, scaleRad2);
  return (psi);//psi is already the negative of the potential
}

/* the profiles return the negative of the result because max_finder()
 *is ACTUALLY meant for finding minima. max_finder() returns the negative 
 * of the function value at the minima found, thus being the max of the function*/
static inline real profile_rho(real r, real * args)
{
  real mass1 = args[0];
  real mass2 = args[1];
  real scaleRad1 = args[2];
  real scaleRad2 = args[3];
  
  real result =  r*r*density( r, mass1, mass2, scaleRad1, scaleRad2);  
  return -result;
}
  
static inline real profile_dist(real v, real * args )
{
  real mass1 = args[0];
  real mass2 = args[1];
  real scaleRad1 = args[2];
  real scaleRad2 = args[3];
  real r= args[4];
  real upperlimit_r=args[5];
  real energy= args[6];
  real result =  v*v*dist_fun(mass1, mass2, scaleRad1, scaleRad2, energy,upperlimit_r);  
  return -result;
}

/*this is a maxfinding routine to find the maximum of the density.
 * It uses Golden Section Search as outlined in Numerical Recipes 3rd edition
 */

static inline real max_finder(real (*profile)(real , real*), real* profileParams, real a, real b, real c, int limit, real tolerance)
{
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

  profile_x1 = (*profile)(x1, profileParams);
  profile_x2 = (*profile)(x2, profileParams);
  
  while (mw_fabs(x3 - x0) > (tolerance * (mw_fabs(x1) + mw_fabs(x2)) ) )
    {
      counter++;
      if (profile_x2 < profile_x1)
	{
	  x0 = x1;
	  x1 = x2;
	  x2 = RATIO * x2 + RATIO_COMPLEMENT * x3;
	  profile_x1 = (real)profile_x2;
	  profile_x2 = (*profile)(x2,profileParams);
	}
      else
	{
	  x3 = x2;
	  x2 = x1;
	  x1 = RATIO * x1 + RATIO_COMPLEMENT * x0;
	  profile_x2 = (real)profile_x1;
	  profile_x1 = (*profile)(x1, profileParams);
	}
	if(counter>limit){break;}
    }

  if (profile_x1 < profile_x2)
    {
      return (-profile_x1);
    }
  else
    {
      return (-profile_x2);
    }
}



/* assigns angles. Allows for non-circular orbits.*/
static inline mwvector angles(dsfmt_t* dsfmtState, real rad)
{
    mwvector vec;
    real phi, theta;

    /*defining some angles*/
    theta = mw_acos( mwXrandom(dsfmtState, -1.0, 1.0) );
    phi =   mwXrandom( dsfmtState, 0.0, 2.0*M_PI );

    /*this is standard formula for x,y,z components in spherical*/
    X(vec) = rad*sin( theta )*cos( phi );    /*x component*/
    Y(vec) = rad*sin( theta )*sin( phi );    /*y component*/
    Z(vec) = rad*cos( theta );               /*z component*/

    return vec;
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
      val = r*r * density(r,  mass1,  mass2,  scaleRad1,  scaleRad2);
  
      if (val/rho_max > u)
      {
       	GOOD_RADIUS = 1;
      }
    }
  return r;
}



static inline real vel_mag(dsfmt_t* dsfmtState,real r, real mass1, real mass2, real scaleRad1, real scaleRad2)
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
  int counter=0;
  real val,v,u,d;
  real energy;
  mwbool returning=FALSE;
  real upperlimit_r;
  
  real mass_en1= mass_en(r, mass1, scaleRad1);
  real mass_en2= mass_en(r, mass2, scaleRad2);
  real v_esc= mw_sqrt( mw_fabs(2.0* (mass_en1+mass_en2)/r));
  
  real parameters[4]= {mass1, mass2, scaleRad1, scaleRad2};
  energy= potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*v_esc*v_esc;
  
  /*want to make sure the root finder returns something proper*/
    do
    {
      upperlimit_r=findRoot(profile_psi,parameters, energy, 0.0, 2.0*(scaleRad1+scaleRad2), dsfmtState); 
      if(isinf(upperlimit_r)==FALSE && upperlimit_r>=0.0 && isnan(upperlimit_r)==FALSE){break;}
      counter++;
      if(counter>=10)
      {
	returning=TRUE;
	break;
      }
    }
    while(1);
    if(returning==TRUE){return 0.0;}
    
    real args[7]={mass1,mass2, scaleRad1, scaleRad2, r, upperlimit_r, energy};
    real dist_max=max_finder(profile_dist, args, 0.0, .5*v_esc, v_esc, 10, 1e-2);
    
    while(1)
    {

      v = (real)mwXrandom(dsfmtState,0.0, v_esc);
      u = (real)mwXrandom(dsfmtState,0.0,1.0);
      
      energy= potential( r, mass1, mass2, scaleRad1, scaleRad2)-0.5*v*v;
      upperlimit_r=findRoot(profile_psi,parameters, energy, 0.0, 2.0*(scaleRad1+scaleRad2), dsfmtState); 
      if(isinf(upperlimit_r)==TRUE || upperlimit_r<=0.0 || isnan(upperlimit_r)==TRUE)
      {continue;}//no point in doing a do-while loophere since we are already in one.
      d=dist_fun(mass1,  mass2,  scaleRad1,  scaleRad2, energy, upperlimit_r);
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
    
    real half_bodies= 0.5*nbody;
    real light_count=0;//counter for the number of light particles assigned
    real mass_light_particle = mass1 / (half_bodies);//half the particles are light matter
    real mass_dark_particle = mass2 / (half_bodies);//half dark matter
    
    /*dark matter type is TRUE or 1. Light matter type is False, or 0*/
    mwbool isdark = TRUE;
    mwbool islight = FALSE;
    int dark= 1;
    int light=0;
    int N=nbody;//integer number of bodies
    real max_light_density;//the max of the light matter density
    real all_r[N];//array to store the radii
    real mass_type[N];//array to store body type
    
    /*the following calculates the massratio if needed:*/
//     real light_needed= (mass1/mass)*nbody;
    
    /*getting the maximum of the density depending on the scale radii*/
    real args[4]= {mass1,mass2, radiusScale1, radiusScale2};
    real rho_max;
    
    
    /*what is happening here: 
     * the r^2*rho function has 2 maxima and 1 minima. We use Newton Raph to 
     * find the max of the density function for each individual component. We then 
     * use the root for those as a starting point for NR to find the max of the combined
     * function and compare to see which is larger. 
     */
    real parameters_light[4]= {mass1, 0.0, radiusScale1, radiusScale2};
    real parameters_dark[4] = {0.0, mass2, radiusScale1, radiusScale2};

    /*finding the max of the individual components*/
    real rho_root_light=findRoot(profile_deriv_rho, parameters_light, 0.0, 0.0,radiusScale1+radiusScale2, prng);
    real rho_max_light= sqr(rho_root_light)*density(rho_root_light,mass1, 0.0, radiusScale1, radiusScale2);
    
    real rho_root_dark=findRoot(profile_deriv_rho, parameters_dark, 0.0, 0.0,radiusScale1 + radiusScale2, prng);
    real rho_max_dark= sqr(rho_root_dark)*density(rho_root_dark, 0.0, mass2, radiusScale1, radiusScale2);
    
    /*finding the max of the combined function starting at previous found roots*/
    real rho_root_lightroot=findRoot(profile_deriv_rho, args, 0.0, rho_max_light,rho_max_light, prng);
    real rho_max_lightroot= sqr(rho_root_lightroot)* density(rho_root_lightroot, mass1, mass2, radiusScale1, radiusScale2);
    
    real rho_root_darkroot=findRoot(profile_deriv_rho, args, 0.0, rho_max_dark,rho_max_dark, prng);
    real rho_max_darkroot= sqr(rho_root_darkroot)* density(rho_root_darkroot, mass1, mass2, radiusScale1, radiusScale2);
    
    
    if(rho_max_lightroot>rho_max_darkroot)
    {rho_max=rho_max_lightroot;}
    else if(rho_max_lightroot<rho_max_darkroot)
    {rho_max=rho_max_darkroot;}
    else
    {rho_max=max_finder(profile_rho, args, 0,radiusScale2, (radiusScale1 + radiusScale2), 20, 1e-4);}
    
    
    memset(&b, 0, sizeof(b));
    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);	
    
      /*getting the radii and velocities for the bodies*/
      for (i = 0; i < nbody; i++)
      {
	  do
	  {
	    r= r_mag(prng, mass1, mass2, radiusScale1, radiusScale2, rho_max);
	    /*to ensure that r is finite and nonzero*/
	    if(isinf(r)==FALSE && r!=0.0 && isnan(r)==FALSE){break;}
	  }
	  while (1);
	  

	  /*storing r's, mass type*/
	  all_r[i]=r;
	  mass_type[i]= dark; //starting them all off dark
      }
      
      /*testing for light matter*/
      i=0;
      light_count=0;
      real coeff= (3.0/2.0)*(minusfivehalves( (3.0/5.0 ) ) );
      real u;
      while(light_count<half_bodies)//only want half the bodies light matter
      {

	if(mass_type[i]==dark)
	{
	  
	  r=all_r[i];
	  /*max_light_density is equal to r^2*rho(r)/ (r^2*rho(r))_max*/
	  max_light_density=coeff* sqr(r)/sqr(radiusScale1)* minusfivehalves( (1.0 + sqr(r)/sqr(radiusScale1)) )  ;
	  u= (real)mwXrandom(prng,0.0,1.0);

	  if( max_light_density < u)
	  {
	    mass_type[i]=light;
	    light_count++;
	  }
	}
	
	if(i>= nbody-1){i=0;}
	else{i++;}
      }
      
      /*this actually gets the position and velocity vectors and pushes table of bodies*/

      for (i = 0; i < nbody; i++)
      {
	  r=all_r[i];
	  
	  if(mass_type[i] == light)
	  {
	    b.bodynode.mass = mass_light_particle;
	    b.bodynode.type = BODY(islight);
	  }
	  else if(mass_type[i]==dark)
	  {
	    b.bodynode.mass = mass_dark_particle;
	    b.bodynode.type = BODY(isdark);
	  }
	  
// 	  mw_printf("\r velocity of particle %i", i);
	  do
	  {
	    v= vel_mag(prng, r, mass1, mass2, radiusScale1, radiusScale2);
	    if(isinf(v)==FALSE && v!=0.0 && isnan(v)==FALSE){break;}
	  }
	  while (1);
	  
  
	  b.vel = vel_vec(prng,  vShift, v);
	  b.bodynode.pos = r_vec(prng, rShift, r);
	  
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


